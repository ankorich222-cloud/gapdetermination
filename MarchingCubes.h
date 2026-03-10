#pragma once
#include "VolumeData.h"
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <cmath>
#include <array>

#include "mc_tables.h"  

class MarchingCubes {
public:
    // 主入口：从体数据提取等值面
    static vtkSmartPointer<vtkPolyData> extract(
        const VolumeData& vol, float isoValue)
    {
        auto points = vtkSmartPointer<vtkPoints>::New();
        auto polys  = vtkSmartPointer<vtkCellArray>::New();
        auto normals= vtkSmartPointer<vtkFloatArray>::New();
        normals->SetNumberOfComponents(3);
        normals->SetName("Normals");

        // 遍历所有 cube（每个 cube 由8个相邻体素构成）
        for (int z = 0; z < vol.dims[2] - 1; ++z)
        for (int y = 0; y < vol.dims[1] - 1; ++y)
        for (int x = 0; x < vol.dims[0] - 1; ++x) {
            processCube(vol, x, y, z, isoValue,
                        points, polys, normals);
        }

        auto polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->SetPolys(polys);
        polyData->GetPointData()->SetNormals(normals);
        return polyData;
    }

private:
    // ── 处理单个 Cube ────────────────────────────────────────────
    static void processCube(
        const VolumeData& vol,
        int x, int y, int z, float iso,
        vtkPoints* pts, vtkCellArray* cells,
        vtkFloatArray* norms)
    {
        // cube 8个顶点的灰度值
        // 顶点编号遵循标准 MC 规范：
        //   4---5
        //  /|  /|
        // 7---6 |
        // | 0-|-1
        // |/  |/
        // 3---2
        float val[8] = {
            vol.at(x,   y,   z  ),  // 0
            vol.at(x+1, y,   z  ),  // 1
            vol.at(x+1, y+1, z  ),  // 2
            vol.at(x,   y+1, z  ),  // 3
            vol.at(x,   y,   z+1),  // 4
            vol.at(x+1, y,   z+1),  // 5
            vol.at(x+1, y+1, z+1),  // 6
            vol.at(x,   y+1, z+1)   // 7
        };

        // 计算 cube 索引（8位，每位=该顶点是否在等值面内侧）
        // 注意：这里假设 val < iso 表示在内侧，这取决于具体应用背景
        int cubeIndex = 0;
        for (int i = 0; i < 8; ++i)
            if (val[i] < iso) cubeIndex |= (1 << i);

        // 全在内/外侧，无交叉
        if (edgeTable[cubeIndex] == 0) return;

        // 12条边上的插值顶点坐标
        // 物理坐标 = 体素坐标 × spacing + origin
        auto voxelToWorld = [&](float fx, float fy, float fz) 
            -> std::array<float, 3> {
            return {
                static_cast<float>(vol.origin[0] + fx * vol.spacing[0]),
                static_cast<float>(vol.origin[1] + fy * vol.spacing[1]),
                static_cast<float>(vol.origin[2] + fz * vol.spacing[2])
            };
        };

        // 12条边端点的体素相对偏移
        static const int edgeVerts[12][2][3] = {
            {{0,0,0},{1,0,0}}, {{1,0,0},{1,1,0}},
            {{0,1,0},{1,1,0}}, {{0,0,0},{0,1,0}},
            {{0,0,1},{1,0,1}}, {{1,0,1},{1,1,1}},
            {{0,1,1},{1,1,1}}, {{0,0,1},{0,1,1}},
            {{0,0,0},{0,0,1}}, {{1,0,0},{1,0,1}},
            {{1,1,0},{1,1,1}}, {{0,1,0},{0,1,1}}
        };

        // 插值计算边上的交点
        std::array<float,3> edgePts[12];
        for (int e = 0; e < 12; ++e) {
            if (!(edgeTable[cubeIndex] & (1 << e))) continue;

            const int* v0_off = edgeVerts[e][0];
            const int* v1_off = edgeVerts[e][1];

            float val0 = vol.at(x + v0_off[0], y + v0_off[1], z + v0_off[2]);
            float val1 = vol.at(x + v1_off[0], y + v1_off[1], z + v1_off[2]);

            // 线性插值：t = (iso - val0) / (val1 - val0)
            float t = (std::abs(val1 - val0) < 1e-6f) ? 0.5f :
                      (iso - val0) / (val1 - val0);

            float fx = x + v0_off[0] + t * (v1_off[0] - v0_off[0]);
            float fy = y + v0_off[1] + t * (v1_off[1] - v0_off[1]);
            float fz = z + v0_off[2] + t * (v1_off[2] - v0_off[2]);

            edgePts[e] = voxelToWorld(fx, fy, fz);
        }

        // 按查找表生成三角形
        for (int i = 0; triTable[cubeIndex][i] != -1; i += 3) {
            int e0 = triTable[cubeIndex][i];
            int e1 = triTable[cubeIndex][i+1];
            int e2 = triTable[cubeIndex][i+2];

            auto& p0 = edgePts[e0];
            auto& p1 = edgePts[e1];
            auto& p2 = edgePts[e2];

            vtkIdType ids[3];
            ids[0] = pts->InsertNextPoint(p0.data());
            ids[1] = pts->InsertNextPoint(p1.data());
            ids[2] = pts->InsertNextPoint(p2.data());

            // 计算面法线
            float v10[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
            float v20[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
            
            float nx = v10[1]*v20[2] - v10[2]*v20[1];
            float ny = v10[2]*v20[0] - v10[0]*v20[2];
            float nz = v10[0]*v20[1] - v10[1]*v20[0];
            
            float len = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (len > 1e-6f) { nx/=len; ny/=len; nz/=len; }

            for (int k = 0; k < 3; ++k)
                norms->InsertNextTuple3(nx, ny, nz);

            cells->InsertNextCell(3, ids);
        }
    }
};
