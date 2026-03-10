#pragma once
#include "VolumeData.h"
#include <vector>
#include <cstdint>
#include <queue>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkMarchingCubes.h>
#include <vtkImageImport.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 增强型并查集
struct UnionFind {
    static int find(int* parent, int i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    }
    static void unite(int* parent, int i, int j) {
        int root_i = find(parent, i);
        int root_j = find(parent, j);
        if (root_i != root_j) parent[root_i] = root_j;
    }
};

struct VoidRegion {
    int id = 0;
    size_t voxelCount = 0;
    double volumeMM3 = 0.0;
    double equivalentDiameterMM = 0.0;
    std::array<double, 3> centroidMM = {0, 0, 0};
    std::array<int, 6> bbox = {0, 0, 0, 0, 0, 0};
};

struct VoidDetectionParams {
    float grayMin = 0.0f;
    float grayMax = 0.0f;
    double minVolumeMM3 = 0.01;
    // ── 核心向量参数 ──
    float angleThresholdDeg = 30.0f; 
    int   tensorWindowSize = 1;      // 张量分析窗口半径 (1 = 3x3x3)
};

class VoidAnalysis {
public:
    // ── 计算结构张量主特征向量 ───────────────
    // 利用窗口内梯度分布求本征方向
    static std::array<float, 3> getPrincipalDirection(const VolumeData& vol, int x, int y, int z, int window) {
        float m11=0, m12=0, m13=0, m22=0, m23=0, m33=0;
        
        for(int dz=-window; dz<=window; ++dz)
        for(int dy=-window; dy<=window; ++dy)
        for(int dx=-window; dx<=window; ++dx) {
            float gx = (vol.at(x+dx+1, y+dy, z+dz) - vol.at(x+dx-1, y+dy, z+dz)) * 0.5f;
            float gy = (vol.at(x+dx, y+dy+1, z+dz) - vol.at(x+dx, y+dy-1, z+dz)) * 0.5f;
            float gz = (vol.at(x+dx, y+dy, z+dz+1) - vol.at(x+dx, y+dy, z+dz-1)) * 0.5f;
            m11 += gx*gx; m12 += gx*gy; m13 += gx*gz;
            m22 += gy*gy; m23 += gy*gz; m33 += gz*gz;
        }
        
        // 简化版特征向量估计：使用最强梯度分量
        // 在工业级实现中，此处应进行 Jacobi 特征值分解
        // 为了实时性，我们计算重心梯度作为方向表示
        float len = std::sqrt(m11 + m22 + m33 + 1e-9f);
        return { std::sqrt(m11)/len, std::sqrt(m22)/len, std::sqrt(m33)/len };
    }

    // ── Step 4 & 5 ──
    static std::vector<uint8_t> createInteriorMaskVoxel(const VolumeData& vol, float isoValue) {
        int dx = vol.dims[0]; int dy = vol.dims[1]; int dz = vol.dims[2];
        size_t total = (size_t)dx * dy * dz;
        std::vector<uint8_t> mask(total, 0);
        std::queue<size_t> q;
        auto add_seed = [&](int x, int y, int z) {
            size_t idx = (size_t)x + (size_t)y * dx + (size_t)z * dx * dy;
            if (vol.at(x, y, z) < isoValue && mask[idx] == 0) { mask[idx] = 1; q.push(idx); }
        };
        for (int z : {0, dz-1}) for (int y = 0; y < dy; ++y) for (int x = 0; x < dx; ++x) add_seed(x, y, z);
        for (int y : {0, dy-1}) for (int z = 0; z < dz; ++z) for (int x = 0; x < dx; ++x) add_seed(x, y, z);
        for (int x : {0, dx-1}) for (int z = 0; z < dz; ++z) for (int y = 0; y < dy; ++y) add_seed(x, y, z);
        const int neighbors[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
        while (!q.empty()) {
            size_t currIdx = q.front(); q.pop();
            int cz = (int)(currIdx / (dx * dy)); int cy = (int)((currIdx / dx) % dy); int cx = (int)(currIdx % dx);
            for (int i = 0; i < 6; ++i) {
                int nx = cx + neighbors[i][0]; int ny = cy + neighbors[i][1]; int nz = cz + neighbors[i][2];
                if (nx >= 0 && nx < dx && ny >= 0 && ny < dy && nz >= 0 && nz < dz) {
                    size_t nIdx = (size_t)nx + (size_t)ny * dx + (size_t)nz * dx * dy;
                    if (mask[nIdx] == 0 && vol.at(nx, ny, nz) < isoValue) { mask[nIdx] = 1; q.push(nIdx); }
                }
            }
        }
        #pragma omp parallel for
        for (long long i = 0; i < (long long)total; ++i) mask[i] = (mask[i] == 1) ? 0 : 1;
        return mask;
    }

    static std::vector<uint8_t> extractCandidateVoids(const VolumeData& vol, const std::vector<uint8_t>& interiorMask, const VoidDetectionParams& params) {
        size_t total = (size_t)vol.dims[0] * vol.dims[1] * vol.dims[2];
        std::vector<uint8_t> candidates(total, 0);
        #pragma omp parallel for
        for (long long i = 0; i < (long long)total; ++i) {
            if (interiorMask[i] > 0 && vol.voxels[i] >= params.grayMin && vol.voxels[i] <= params.grayMax) candidates[i] = 1;
        }
        return candidates;
    }

    // ── 向量计算驱动的连通域分析  ──────────────────
    static std::vector<VoidRegion> labelAndAnalyze(
        const VolumeData& vol,
        std::vector<uint8_t>& mask,
        const VoidDetectionParams& params,
        std::vector<int>& outLabelVol) 
    {
        int dx = vol.dims[0]; int dy = vol.dims[1]; int dz = vol.dims[2];
        size_t total = (size_t)dx * dy * dz;
        outLabelVol.assign(total, -1);

        float cosThreshold = std::cos(params.angleThresholdDeg * (float)M_PI / 180.0f);

        // 1. 预计算向量场 
        printf("Step A: Computing structure tensor field...\n");
        std::vector<std::array<float, 3>> directions(total);
        #pragma omp parallel for
        for (long long i = 0; i < (long long)total; ++i) {
            if (mask[i] == 0) continue;
            int cz = (int)(i / (dx * dy)); int cy = (int)((i / dx) % dy); int cx = (int)(i % dx);
            directions[i] = getPrincipalDirection(vol, cx, cy, cz, params.tensorWindowSize);
        }

        // 2. 场驱动合并 (并查集扫描线)
        printf("Step B: Directional field merging...\n");
        const int nOffsets[13][3] = {
            {-1,0,0}, {0,-1,0}, {0,0,-1}, {-1,-1,0}, {-1,1,0}, {-1,0,-1}, {-1,0,1}, {0,-1,-1}, {0,-1,1},
            {-1,-1,-1}, {-1,-1,1}, {-1,1,-1}, {-1,1,1}
        };

        for (int z = 0; z < dz; ++z) {
            for (int y = 0; y < dy; ++y) {
                for (int x = 0; x < dx; ++x) {
                    size_t curr = (size_t)x + (size_t)y * dx + (size_t)z * dx * dy;
                    if (mask[curr] == 0) continue;
                    outLabelVol[curr] = (int)curr;

                    auto& dirCurr = directions[curr];
                    for (int i = 0; i < 13; ++i) {
                        int nx = x + nOffsets[i][0]; int ny = y + nOffsets[i][1]; int nz = z + nOffsets[i][2];
                        if (nx >= 0 && nx < dx && ny >= 0 && ny < dy && nz >= 0 && nz < dz) {
                            size_t nIdx = (size_t)nx + (size_t)ny * dx + (size_t)nz * dx * dy;
                            if (outLabelVol[nIdx] != -1) {
                                auto& dirNext = directions[nIdx];
                                float dot = std::abs(dirCurr[0]*dirNext[0] + dirCurr[1]*dirNext[1] + dirCurr[2]*dirNext[2]);
                                if (dot >= cosThreshold) UnionFind::unite(outLabelVol.data(), (int)curr, (int)nIdx);
                            }
                        }
                    }
                }
            }
        }

        // 3. ID 扁平化与统计
        printf("Step C: Finalizing results...\n");
        std::map<int, int> rootMap;
        std::vector<VoidRegion> regions;
        int nextId = 1;
        for (size_t i = 0; i < total; ++i) {
            if (outLabelVol[i] == -1) { outLabelVol[i] = 0; continue; }
            int root = UnionFind::find(outLabelVol.data(), (int)i);
            if (rootMap.find(root) == rootMap.end()) {
                rootMap[root] = nextId++;
                VoidRegion r; r.id = rootMap[root]; regions.push_back(r);
            }
            int fId = rootMap[root];
            outLabelVol[i] = fId;
            regions[fId-1].voxelCount++;
        }

        double vVol = vol.spacing[0] * vol.spacing[1] * vol.spacing[2];
        std::vector<VoidRegion> filtered;
        for (auto& r : regions) {
            r.volumeMM3 = r.voxelCount * vVol;
            if (r.volumeMM3 >= params.minVolumeMM3) {
                r.equivalentDiameterMM = std::pow((6.0 * r.volumeMM3) / M_PI, 1.0 / 3.0);
                filtered.push_back(r);
            }
        }
        return filtered;
    }

    static bool saveResults(const std::string& baseName, const std::vector<VoidRegion>& voids, const std::vector<int>& labelVol, const VolumeData& vol) {
        std::string csvPath = baseName + "_voids.csv";
        FILE* fp = fopen(csvPath.c_str(), "w");
        if (fp) {
            fprintf(fp, "ID,VoxelCount,Volume(mm3),EquivDiameter(mm)\n");
            for (const auto& v : voids) fprintf(fp, "%d,%llu,%.4f,%.4f\n", v.id, (unsigned long long)v.voxelCount, v.volumeMM3, v.equivalentDiameterMM);
            fclose(fp);
        }
        std::string rawPath = baseName + "_label.raw";
        FILE* fr = fopen(rawPath.c_str(), "wb");
        if (fr) { fwrite(labelVol.data(), sizeof(int), labelVol.size(), fr); fclose(fr); }
        return true;
    }

    static vtkSmartPointer<vtkPolyData> generateVoidMesh(const std::vector<int>& labelVol, const VolumeData& vol) {
        auto importer = vtkSmartPointer<vtkImageImport>::New();
        importer->SetDataScalarTypeToInt(); importer->SetNumberOfScalarComponents(1);
        importer->SetWholeExtent(0, vol.dims[0]-1, 0, vol.dims[1]-1, 0, vol.dims[2]-1);
        importer->SetDataExtentToWholeExtent(); importer->SetDataSpacing(vol.spacing[0], vol.spacing[1], vol.spacing[2]);
        importer->SetDataOrigin(vol.origin[0], vol.origin[1], vol.origin[2]);
        importer->SetImportVoidPointer(const_cast<int*>(labelVol.data())); importer->Update();
        auto mc = vtkSmartPointer<vtkMarchingCubes>::New();
        mc->SetInputConnection(importer->GetOutputPort()); mc->SetValue(0, 0.5); mc->Update();
        return mc->GetOutput();
    }
};
