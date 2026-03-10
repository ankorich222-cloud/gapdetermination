#pragma once
#include "VolumeData.h"
#include <vector>
#include <cstdint>
#include <queue>
#include <algorithm>
#include <array>
#include <cmath>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkMarchingCubes.h>
#include <vtkImageImport.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 单个空隙的几何描述
struct VoidRegion {
    int id = 0;
    size_t voxelCount = 0;
    double volumeMM3 = 0.0;
    double equivalentDiameterMM = 0.0;
    std::array<double, 3> centroidMM = {0, 0, 0};
    std::array<int, 6> bbox = {0, 0, 0, 0, 0, 0}; // [minX, maxX, minY, maxY, minZ, maxZ]
    std::array<int, 3> seedVoxel = {0, 0, 0};    // 种子点坐标
};

struct VoidDetectionParams {
    float grayMin = 0.0f;
    float grayMax = 0.0f;
    int   connectivity = 26; 
    double minVolumeMM3 = 0.01; // 过滤掉过小的噪声点
};

class VoidAnalysis {
public:
    // ── Step 4: 提取内部掩膜 ─────────────────────────────────────
    static std::vector<uint8_t> createInteriorMaskVoxel(const VolumeData& vol, float isoValue) {
        int dx = vol.dims[0]; int dy = vol.dims[1]; int dz = vol.dims[2];
        size_t total = (size_t)dx * dy * dz;
        std::vector<uint8_t> mask(total, 0);
        std::queue<size_t> q;

        auto add_seed = [&](int x, int y, int z) {
            size_t idx = (size_t)x + (size_t)y * dx + (size_t)z * dx * dy;
            if (vol.at(x, y, z) < isoValue && mask[idx] == 0) {
                mask[idx] = 1; q.push(idx);
            }
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

    // ── Step 5: 提取候选空隙 ─────────────────────────────────────
    static std::vector<uint8_t> extractCandidateVoids(const VolumeData& vol, const std::vector<uint8_t>& interiorMask, const VoidDetectionParams& params) {
        size_t total = (size_t)vol.dims[0] * vol.dims[1] * vol.dims[2];
        std::vector<uint8_t> candidates(total, 0);
        #pragma omp parallel for
        for (long long i = 0; i < (long long)total; ++i) {
            if (interiorMask[i] > 0 && vol.voxels[i] >= params.grayMin && vol.voxels[i] <= params.grayMax) candidates[i] = 1;
        }
        return candidates;
    }

    // ── Step 6 & 7: 连通域分割与空隙定量分析 ──────────────────────
    static std::vector<VoidRegion> labelAndAnalyze(
        const VolumeData& vol,
        std::vector<uint8_t>& candidateMask,
        const VoidDetectionParams& params,
        std::vector<int>& outLabelVol) 
    {
        int dx = vol.dims[0]; int dy = vol.dims[1]; int dz = vol.dims[2];
        size_t total = (size_t)dx * dy * dz;
        double voxelVol = vol.spacing[0] * vol.spacing[1] * vol.spacing[2];

        outLabelVol.assign(total, 0);
        std::vector<VoidRegion> regions;
        std::vector<uint8_t> visited(total, 0);
        int nextId = 1;

        std::vector<std::array<int, 3>> neighbors;
        for (int sz = -1; sz <= 1; ++sz)
        for (int sy = -1; sy <= 1; ++sy)
        for (int sx = -1; sx <= 1; ++sx) {
            if (sx == 0 && sy == 0 && sz == 0) continue;
            neighbors.push_back({sx, sy, sz});
        }

        for (int z = 0; z < dz; ++z)
        for (int y = 0; y < dy; ++y)
        for (int x = 0; x < dx; ++x) {
            size_t idx = (size_t)x + (size_t)y * dx + (size_t)z * dx * dy;
            if (candidateMask[idx] == 1 && !visited[idx]) {
                VoidRegion region;
                region.id = nextId++;
                region.seedVoxel = {x, y, z};
                region.bbox = {x, x, y, y, z, z};
                double sumX = 0, sumY = 0, sumZ = 0;
                std::vector<size_t> currentRegionVoxels;

                std::queue<size_t> q;
                q.push(idx);
                visited[idx] = 1;

                while (!q.empty()) {
                    size_t curr = q.front(); q.pop();
                    currentRegionVoxels.push_back(curr);
                    int cz = (int)(curr / (dx * dy)); int cy = (int)((curr / dx) % dy); int cx = (int)(curr % dx);
                    region.voxelCount++;
                    sumX += cx; sumY += cy; sumZ += cz;
                    region.bbox[0] = std::min(region.bbox[0], cx); region.bbox[1] = std::max(region.bbox[1], cx);
                    region.bbox[2] = std::min(region.bbox[2], cy); region.bbox[3] = std::max(region.bbox[3], cy);
                    region.bbox[4] = std::min(region.bbox[4], cz); region.bbox[5] = std::max(region.bbox[5], cz);

                    for (auto& n : neighbors) {
                        int nx = cx + n[0]; int ny = cy + n[1]; int nz = cz + n[2];
                        if (nx >= 0 && nx < dx && ny >= 0 && ny < dy && nz >= 0 && nz < dz) {
                            size_t nIdx = (size_t)nx + (size_t)ny * dx + (size_t)nz * dx * dy;
                            if (candidateMask[nIdx] == 1 && !visited[nIdx]) {
                                visited[nIdx] = 1; q.push(nIdx);
                            }
                        }
                    }
                }

                region.volumeMM3 = region.voxelCount * voxelVol;
                if (region.volumeMM3 >= params.minVolumeMM3) {
                    for (size_t vIdx : currentRegionVoxels) outLabelVol[vIdx] = region.id;
                    region.centroidMM[0] = (sumX / region.voxelCount) * vol.spacing[0] + vol.origin[0];
                    region.centroidMM[1] = (sumY / region.voxelCount) * vol.spacing[1] + vol.origin[1];
                    region.centroidMM[2] = (sumZ / region.voxelCount) * vol.spacing[2] + vol.origin[2];
                    region.equivalentDiameterMM = std::pow((6.0 * region.volumeMM3) / M_PI, 1.0/3.0);
                    regions.push_back(region);
                } else {
                    nextId--; 
                }
            }
        }
        return regions;
    }

    // ── Step 10: 结果导出 ─────────────────────────────────────────
    static bool saveResults(
        const std::string& baseName,
        const std::vector<VoidRegion>& voids,
        const std::vector<int>& labelVol,
        const VolumeData& vol)
    {
        std::string csvPath = baseName + "_voids.csv";
        FILE* fp = fopen(csvPath.c_str(), "w");
        if (fp) {
            fprintf(fp, "ID,VoxelCount,Volume(mm3),EquivDiameter(mm),CentroidX,CentroidY,CentroidZ,minX,maxX,minY,maxY,minZ,maxZ\n");
            for (const auto& v : voids) {
                fprintf(fp, "%d,%llu,%.4f,%.4f,%.3f,%.3f,%.3f,%d,%d,%d,%d,%d,%d\n",
                    v.id, (unsigned long long)v.voxelCount, v.volumeMM3, v.equivalentDiameterMM,
                    v.centroidMM[0], v.centroidMM[1], v.centroidMM[2],
                    v.bbox[0], v.bbox[1], v.bbox[2], v.bbox[3], v.bbox[4], v.bbox[5]);
            }
            fclose(fp);
            printf("Saved: %s\n", csvPath.c_str());
        }

        std::string rawPath = baseName + "_label.raw";
        FILE* fr = fopen(rawPath.c_str(), "wb");
        if (fr) {
            fwrite(labelVol.data(), sizeof(int), labelVol.size(), fr);
            fclose(fr);
            printf("Saved: %s (Format: Int32, Dims: %d x %d x %d)\n", 
                   rawPath.c_str(), vol.dims[0], vol.dims[1], vol.dims[2]);
        }
        return true;
    }

    // ── Step 9.3: 将空隙标签体转换为 3D 网格 ───────────────────────
    static vtkSmartPointer<vtkPolyData> generateVoidMesh(
        const std::vector<int>& labelVol,
        const VolumeData& vol)
    {
        auto importer = vtkSmartPointer<vtkImageImport>::New();
        importer->SetDataScalarTypeToInt();
        importer->SetNumberOfScalarComponents(1);
        importer->SetWholeExtent(0, vol.dims[0]-1, 0, vol.dims[1]-1, 0, vol.dims[2]-1);
        importer->SetDataExtentToWholeExtent();
        importer->SetDataSpacing(vol.spacing[0], vol.spacing[1], vol.spacing[2]);
        importer->SetDataOrigin(vol.origin[0], vol.origin[1], vol.origin[2]);
        importer->SetImportVoidPointer(const_cast<int*>(labelVol.data()));
        importer->Update();

        auto mc = vtkSmartPointer<vtkMarchingCubes>::New();
        mc->SetInputConnection(importer->GetOutputPort());
        mc->SetValue(0, 0.5); // 区分 0 (背景) 和 >=1 (空隙)
        mc->Update();

        return mc->GetOutput();
    }
};
