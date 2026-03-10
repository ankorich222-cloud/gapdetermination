#pragma once
#include "VolumeData.h"
#include <opencv2/opencv.hpp>
#include <algorithm>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkSmoothPolyDataFilter.h>

struct SurfaceParams {
    float background;   // 背景灰度值上限
    float material;     // 材料灰度值下限
    float isoValue;     // 等值面阈值（表面）
};

// 高级模式参数：法向搜索距离优化边界
struct AdvancedSurfaceParams {
    bool  enabled            = false;
    float normalSearchDistance = 0.0f;   // 搜索距离（单位由 useMillimeter 决定）
    bool  useMillimeter      = false;    // true: mm；false: voxel 近似
    float searchStep         = 0.5f;     // 采样步长（同上单位）
    float maxVertexShift     = 2.0f;     // 顶点最大位移（同上单位）
    float gradientThreshold  = 0.0f;     // 最小梯度阈值（0 表示不限制）
    int   normalSmoothIterations = 0;    // 可选：后续几何平滑迭代次数
};

class SurfaceDetermination {
public:
    // ── 从直方图自动估算背景/材料双峰 ──────────────────────────
    static SurfaceParams autoDetect(const VolumeData& vol,
                                     int numBins = 512) {
        std::vector<float> hist = computeHistogram(vol, numBins);
        float range = vol.maxVal - vol.minVal;

        // 2. 高斯平滑直方图（抑制噪声峰）
        cv::Mat histMat(1, numBins, CV_32F);
        for (int i = 0; i < numBins; ++i)
            histMat.at<float>(i) = hist[i];
        cv::GaussianBlur(histMat, histMat, cv::Size(9, 1), 2.0);

        // 3. 找两个主峰（背景峰 + 材料峰）
        int peak1 = 0, peak2 = 0;
        float val1 = 0, val2 = 0;
        for (int i = 0; i < numBins; ++i) {
            float v = histMat.at<float>(i);
            if (v > val1) { val1 = v; peak1 = i; }
        }
        // 第二峰：排除第一峰附近 10% 范围
        int exclusion = numBins / 10;
        for (int i = 0; i < numBins; ++i) {
            if (std::abs(i - peak1) < exclusion) continue;
            float v = histMat.at<float>(i);
            if (v > val2) { val2 = v; peak2 = i; }
        }
        if (peak2 < peak1) std::swap(peak1, peak2); // 确保 peak1=背景

        // 4. 转换回灰度值
        auto binToVal = [&](int bin) {
            return vol.minVal + bin / float(numBins - 1) * range;
        };

        SurfaceParams p;
        p.background = binToVal(peak1);
        p.material   = binToVal(peak2);
        p.isoValue   = (p.background + p.material) * 0.5f; // ISO-50
        return p;
    }

    // ── 预计算直方图 ─────────────────────────────────────────────
    static std::vector<float> computeHistogram(const VolumeData& vol, int numBins = 512) {
        std::vector<float> hist(numBins, 0.f);
        float range = vol.maxVal - vol.minVal;
        if (range <= 0) return hist;

        for (float v : vol.voxels) {
            int bin = static_cast<int>((v - vol.minVal) / range * (numBins - 1));
            hist[std::clamp(bin, 0, numBins - 1)]++;
        }
        return hist;
    }

    // ── 可视化直方图（使用预计算数据） ─────────────────────────────
    static cv::Mat visualizeHistogram(const VolumeData& vol,
                                       const std::vector<float>& precomputedHist,
                                       const SurfaceParams& params,
                                       int width = 800,
                                       int height = 300) {
        int numBins = (int)precomputedHist.size();
        float maxH = *std::max_element(precomputedHist.begin(), precomputedHist.end());
        float range = vol.maxVal - vol.minVal;

        cv::Mat img(height, width, CV_8UC3, cv::Scalar(30, 30, 30));
        int barW = std::max(1, width / numBins);

        // 绘制直方图柱
        for (int i = 0; i < numBins; ++i) {
            int barH = static_cast<int>(precomputedHist[i] / maxH * (height - 20));
            int x = i * width / numBins;
            cv::rectangle(img,
                cv::Point(x, height - barH),
                cv::Point(x + barW, height),
                cv::Scalar(100, 180, 100), -1);
        }

        // 绘制三条标记线
        auto drawLine = [&](float val, cv::Scalar color, const std::string& label) {
            int x = static_cast<int>((val - vol.minVal) / range * width);
            cv::line(img, {x, 0}, {x, height}, color, 2);
            cv::putText(img, label, {x + 4, 20},
                        cv::FONT_HERSHEY_SIMPLEX, 0.5, color, 1);
        };
        drawLine(params.background, {255, 100, 100}, "BG");
        drawLine(params.material,   {100, 100, 255}, "MAT");
        drawLine(params.isoValue,   {255, 255,   0}, "ISO");

        return img;
    }

    // ── 切片轮廓预览（三方向实时） ───────────────────────────────
    // axis: 0=XY, 1=XZ, 2=YZ
    static cv::Mat sliceWithContour(const VolumeData& vol,
                                     int axis, int sliceIdx,
                                     float isoValue) {
        int w, h;
        if (axis == 0) { w = vol.dims[0]; h = vol.dims[1]; }
        else if (axis == 1) { w = vol.dims[0]; h = vol.dims[2]; }
        else               { w = vol.dims[1]; h = vol.dims[2]; }

        cv::Mat gray(h, w, CV_32F);
        for (int row = 0; row < h; ++row)
        for (int col = 0; col < w; ++col) {
            float v = 0;
            if (axis == 0) v = vol.at(col, row, sliceIdx);
            else if (axis == 1) v = vol.at(col, sliceIdx, row);
            else                v = vol.at(sliceIdx, col, row);
            gray.at<float>(row, col) = v;
        }

        // 归一化到 0-255 显示
        cv::Mat display;
        cv::normalize(gray, display, 0, 255, cv::NORM_MINMAX);
        display.convertTo(display, CV_8U);
        cv::cvtColor(display, display, cv::COLOR_GRAY2BGR);

        // 提取并叠加等值线（轮廓）
        cv::Mat binary;
        cv::threshold(gray, binary, isoValue, 255, cv::THRESH_BINARY);
        binary.convertTo(binary, CV_8U);
        std::vector<std::vector<cv::Point>> contours;
        cv::findContours(binary, contours, cv::RETR_EXTERNAL,
                         cv::CHAIN_APPROX_SIMPLE);
        cv::drawContours(display, contours, -1,
                         cv::Scalar(0, 255, 255), 1); // 黄色轮廓

        return display;
    }

    // ── 高级模式：沿法向搜索优化边界位置 ───────────────────────────
    static vtkSmartPointer<vtkPolyData> refineSurfaceAlongNormals(
        const VolumeData& vol,
        vtkSmartPointer<vtkPolyData> surface,
        const SurfaceParams& surfParams,
        const AdvancedSurfaceParams& adv);

    // ── 2D 预览增强：优化 2D 轮廓点 ──────────────────────────────
    // 用于在切片预览时实时呈现高级模式效果，无需重构 3D Mesh
    static void refineContour2D(
        const VolumeData& vol,
        int axis, int sliceIdx,
        std::vector<std::vector<cv::Point>>& contours,
        const SurfaceParams& surfParams,
        const AdvancedSurfaceParams& adv);

private:
    // 核心算法内核：在给定采样序列中寻找最优边界偏移
    static bool findBestOffsetAlongNormal(
        const std::vector<float>& vals,
        const std::vector<double>& ts,
        float iso,
        float gradThresh,
        double& outBestT) 
    {
        int bestIdx = -1;
        double bestGrad = static_cast<double>(gradThresh);
        size_t N = ts.size();
        for (size_t i = 1; i + 1 < N; ++i) {
            float v0 = vals[i - 1] - iso;
            float v1 = vals[i] - iso;
            float v2 = vals[i + 1] - iso;

            // 检查灰度穿越
            bool hasCross = (v0 * v1 <= 0.f) || (v1 * v2 <= 0.f);
            if (!hasCross) continue;

            double dt = ts[i + 1] - ts[i - 1];
            if (std::abs(dt) < 1e-6) continue;

            // 中心差分求梯度
            double grad = std::abs(static_cast<double>(vals[i + 1] - vals[i - 1]) / dt);
            if (grad > bestGrad) {
                bestGrad = grad;
                bestIdx = static_cast<int>(i);
            }
        }
        if (bestIdx >= 0) {
            outBestT = ts[bestIdx];
            return true;
        }
        return false;
    }
};

// 实现部分
inline vtkSmartPointer<vtkPolyData> SurfaceDetermination::refineSurfaceAlongNormals(
    const VolumeData& vol,
    vtkSmartPointer<vtkPolyData> surface,
    const SurfaceParams& surfParams,
    const AdvancedSurfaceParams& adv)
{
    if (!adv.enabled || !surface) return surface;

    vtkPoints* points = surface->GetPoints();
    vtkDataArray* nArray = surface->GetPointData()->GetNormals();
    if (!points || !nArray) return surface;

    vtkIdType numPts = points->GetNumberOfPoints();
    if (numPts == 0) return surface;

    auto unitToWorld = [&](float v) -> double {
        if (adv.useMillimeter) return static_cast<double>(v);
        double avgSpacing = (vol.spacing[0] + vol.spacing[1] + vol.spacing[2]) / 3.0;
        return static_cast<double>(v) * avgSpacing;
    };

    const double searchDistWorld = unitToWorld(adv.normalSearchDistance);
    const double searchStepWorld = std::max(1e-3, unitToWorld(adv.searchStep));
    const double maxShiftWorld   = unitToWorld(adv.maxVertexShift);
    const float  iso             = surfParams.isoValue;
    const float  gradThresh      = adv.gradientThreshold;

    // 预估采样点数量
    const size_t numSamples = static_cast<size_t>(searchDistWorld * 2 / searchStepWorld) + 5;

    // ── 并行化顶点处理 ──────────────────────────────────────────
    #pragma omp parallel
    {
        // 每个线程私有的临时容器，复用内存并避免竞争
        std::vector<double> ts_private;
        std::vector<float>  vals_private;
        ts_private.reserve(numSamples);
        vals_private.reserve(numSamples);

        #pragma omp for schedule(guided)
        for (vtkIdType pid = 0; pid < numPts; ++pid) {
            double p[3], n[3];
            points->GetPoint(pid, p);
            nArray->GetTuple(pid, n);
            
            double nLen = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
            if (nLen < 1e-6) continue;
            n[0] /= nLen; n[1] /= nLen; n[2] /= nLen;

            ts_private.clear(); 
            vals_private.clear();

            for (double t = -searchDistWorld; t <= searchDistWorld + 1e-7; t += searchStepWorld) {
                ts_private.push_back(t);
                vals_private.push_back(vol.sampleTrilinearWorld(p[0] + t * n[0], p[1] + t * n[1], p[2] + t * n[2]));
            }

            double tBest = 0;
            if (findBestOffsetAlongNormal(vals_private, ts_private, iso, gradThresh, tBest)) {
                if (std::abs(tBest) > maxShiftWorld) tBest = (tBest > 0 ? maxShiftWorld : -maxShiftWorld);
                double newP[3] = { p[0] + tBest * n[0], p[1] + tBest * n[1], p[2] + tBest * n[2] };
                points->SetPoint(pid, newP);
            }
        }
    }

    points->Modified();
    if (adv.normalSmoothIterations > 0) {
        auto smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
        smoother->SetInputData(surface);
        smoother->SetNumberOfIterations(adv.normalSmoothIterations);
        smoother->SetRelaxationFactor(0.1);
        smoother->Update();
        surface = smoother->GetOutput();
    }
    return surface;
}

inline void SurfaceDetermination::refineContour2D(
    const VolumeData& vol,
    int axis, int sliceIdx,
    std::vector<std::vector<cv::Point>>& contours,
    const SurfaceParams& surfParams,
    const AdvancedSurfaceParams& adv)
{
    if (!adv.enabled) return;

    auto unitToVoxel = [&](float v) -> float {
        if (!adv.useMillimeter) return v;
        // 使用该平面内平均 spacing 作为单位换算
        double s0, s1;
        if (axis == 0) { s0 = vol.spacing[0]; s1 = vol.spacing[1]; }
        else if (axis == 1) { s0 = vol.spacing[0]; s1 = vol.spacing[2]; }
        else { s0 = vol.spacing[1]; s1 = vol.spacing[2]; }
        return v / static_cast<float>((s0 + s1) * 0.5);
    };

    float searchDistVoxel = unitToVoxel(adv.normalSearchDistance);
    float searchStepVoxel = std::max(0.1f, unitToVoxel(adv.searchStep));
    float maxShiftVoxel   = unitToVoxel(adv.maxVertexShift);

    std::vector<double> ts;
    std::vector<float> vals;
    ts.reserve(static_cast<size_t>(searchDistVoxel * 2 / searchStepVoxel) + 5);
    vals.reserve(ts.capacity());

    for (auto& poly : contours) {
        for (auto& pt : poly) {
            // 在 2D 切片上计算局部梯度作为搜索方向
            float gx, gy;
            if (axis == 0) { // XY
                gx = (vol.at(pt.x + 1, pt.y, sliceIdx) - vol.at(pt.x - 1, pt.y, sliceIdx)) * 0.5f;
                gy = (vol.at(pt.x, pt.y + 1, sliceIdx) - vol.at(pt.x, pt.y - 1, sliceIdx)) * 0.5f;
            } else if (axis == 1) { // XZ
                gx = (vol.at(pt.x + 1, sliceIdx, pt.y) - vol.at(pt.x - 1, sliceIdx, pt.y)) * 0.5f;
                gy = (vol.at(pt.x, sliceIdx, pt.y + 1) - vol.at(pt.x, sliceIdx, pt.y - 1)) * 0.5f;
            } else { // YZ
                gx = (vol.at(sliceIdx, pt.x + 1, pt.y) - vol.at(sliceIdx, pt.x - 1, pt.y)) * 0.5f;
                gy = (vol.at(sliceIdx, pt.x, pt.y + 1) - vol.at(sliceIdx, pt.x, pt.y - 1)) * 0.5f;
            }

            float len = std::sqrt(gx * gx + gy * gy);
            if (len < 1e-4f) continue;
            float nx = gx / len; float ny = gy / len;

            ts.clear(); vals.clear();
            for (float t = -searchDistVoxel; t <= searchDistVoxel + 1e-4f; t += searchStepVoxel) {
                ts.push_back(static_cast<double>(t));
                float vx = pt.x + t * nx;
                float vy = pt.y + t * ny;
                if (axis == 0)      vals.push_back(vol.sampleTrilinearIndex(vx, vy, (float)sliceIdx));
                else if (axis == 1) vals.push_back(vol.sampleTrilinearIndex(vx, (float)sliceIdx, vy));
                else                vals.push_back(vol.sampleTrilinearIndex((float)sliceIdx, vx, vy));
            }

            double tBest = 0;
            if (findBestOffsetAlongNormal(vals, ts, surfParams.isoValue, adv.gradientThreshold, tBest)) {
                if (std::abs(tBest) > maxShiftVoxel) tBest = (tBest > 0 ? maxShiftVoxel : -maxShiftVoxel);
                pt.x += static_cast<int>(std::round(tBest * nx));
                pt.y += static_cast<int>(std::round(tBest * ny));
            }
        }
    }
}
