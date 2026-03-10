#pragma once
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkDICOMImageReader.h>

class VolumeData {
public:
    // 体素数据（线性存储，顺序 x + y*dimX + z*dimX*dimY）
    std::vector<float> voxels;
    
    std::array<int, 3>    dims;     // [dimX, dimY, dimZ]
    std::array<double, 3> spacing;  // 体素尺寸(mm)
    std::array<double, 3> origin;   // 原点坐标
    
    float minVal, maxVal;           // 灰度范围

    // 从 VTK ImageData 加载
    bool loadFromVTK(vtkSmartPointer<vtkImageData> imageData) {
        imageData->GetDimensions(dims.data());
        imageData->GetSpacing(spacing.data());
        imageData->GetOrigin(origin.data());

        int total = dims[0] * dims[1] * dims[2];
        voxels.resize(total);

        // 统一转换为 float
        for (int i = 0; i < total; ++i) {
            voxels[i] = static_cast<float>(
                imageData->GetScalarComponentAsFloat(
                    i % dims[0],
                    (i / dims[0]) % dims[1],
                    i / (dims[0] * dims[1]), 0));
        }

        auto [minIt, maxIt] = std::minmax_element(voxels.begin(), voxels.end());
        minVal = *minIt;
        maxVal = *maxIt;
        return true;
    }

    // 读取 .mhd 文件
    bool loadMHD(const std::string& path) {
        auto reader = vtkSmartPointer<vtkMetaImageReader>::New();
        reader->SetFileName(path.c_str());
        reader->Update();
        return loadFromVTK(reader->GetOutput());
    }

    // 读取 .raw 文件
    bool loadRAW(const std::string& path, 
                 int dx, int dy, int dz, 
                 double sx, double sy, double sz,
                 const std::string& dataType = "short") 
    {
        FILE* fp = fopen(path.c_str(), "rb");
        if (!fp) return false;

        dims = {dx, dy, dz};
        spacing = {sx, sy, sz};
        origin = {0, 0, 0};

        size_t total = (size_t)dx * dy * dz;
        voxels.resize(total);

        if (dataType == "short" || dataType == "int16") {
            std::vector<short> buffer(total);
            fread(buffer.data(), sizeof(short), total, fp);
            for (size_t i = 0; i < total; ++i) voxels[i] = static_cast<float>(buffer[i]);
        } else if (dataType == "float") {
            fread(voxels.data(), sizeof(float), total, fp);
        } else if (dataType == "uchar" || dataType == "uint8") {
            std::vector<unsigned char> buffer(total);
            fread(buffer.data(), sizeof(unsigned char), total, fp);
            for (size_t i = 0; i < total; ++i) voxels[i] = static_cast<float>(buffer[i]);
        } else {
            fclose(fp);
            return false;
        }

        fclose(fp);
        auto [minIt, maxIt] = std::minmax_element(voxels.begin(), voxels.end());
        minVal = *minIt;
        maxVal = *maxIt;
        return true;
    }

    // 安全取值（含边界检查）
    inline float at(int x, int y, int z) const {
        if (x < 0 || x >= dims[0] ||
            y < 0 || y >= dims[1] ||
            z < 0 || z >= dims[2]) return 0.f;

        //强制使用size_t类型 避免32位溢出
        return voxels[(size_t)x + (size_t)y * dims[0] + (size_t)z * dims[0] * dims[1]];
    }

    // 三线性插值：体素索引坐标空间 (浮点 i,j,k)
    inline float sampleTrilinearIndex(float fx, float fy, float fz) const {
        int x0 = static_cast<int>(std::floor(fx));
        int y0 = static_cast<int>(std::floor(fy));
        int z0 = static_cast<int>(std::floor(fz));
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        int z1 = z0 + 1;

        if (x0 < 0 || y0 < 0 || z0 < 0 ||
            x1 >= dims[0] || y1 >= dims[1] || z1 >= dims[2]) {
            return 0.f;
        }

        float tx = fx - static_cast<float>(x0);
        float ty = fy - static_cast<float>(y0);
        float tz = fz - static_cast<float>(z0);

        float c000 = at(x0, y0, z0);
        float c100 = at(x1, y0, z0);
        float c010 = at(x0, y1, z0);
        float c110 = at(x1, y1, z0);
        float c001 = at(x0, y0, z1);
        float c101 = at(x1, y0, z1);
        float c011 = at(x0, y1, z1);
        float c111 = at(x1, y1, z1);

        float c00 = c000 * (1 - tx) + c100 * tx;
        float c01 = c001 * (1 - tx) + c101 * tx;
        float c10 = c010 * (1 - tx) + c110 * tx;
        float c11 = c011 * (1 - tx) + c111 * tx;

        float c0 = c00 * (1 - ty) + c10 * ty;
        float c1 = c01 * (1 - ty) + c11 * ty;

        return c0 * (1 - tz) + c1 * tz;
    }

    // 三线性插值：世界坐标 (mm) 空间
    inline float sampleTrilinearWorld(double x, double y, double z) const {
        double fx = (x - origin[0]) / spacing[0];
        double fy = (y - origin[1]) / spacing[1];
        double fz = (z - origin[2]) / spacing[2];
        return sampleTrilinearIndex(static_cast<float>(fx),
                                    static_cast<float>(fy),
                                    static_cast<float>(fz));
    }
};
