#include "VolumeData.h"
#include "SurfaceDetermination.h"
#include "MarchingCubes.h"
#include "VoidAnalysis.h"

#include <vtkSmartPointer.h>
#include <vtkSTLWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <opencv2/highgui.hpp>
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.mhd/input.raw>" << std::endl;
        return 1;
    }

    std::string filePath = argv[1];
    VolumeData vol;

    // ── 1. 加载数据 ──────────────────────────────────────────────
    if (filePath.substr(filePath.find_last_of(".") + 1) == "raw") {
        int dx, dy, dz; double sx, sy, sz; std::string dataType;
        std::cout << "Enter dims (X Y Z): "; std::cin >> dx >> dy >> dz;
        std::cout << "Enter spacing: "; std::cin >> sx >> sy >> sz;
        std::cout << "Type (short/float/uchar): "; std::cin >> dataType;
        if (!vol.loadRAW(filePath, dx, dy, dz, sx, sy, sz, dataType)) return 1;
    } else {
        if (!vol.loadMHD(filePath)) return 1;
    }

    // ── 2. 初始化参数 ────────────────────────────────────────────
    SurfaceParams params = SurfaceDetermination::autoDetect(vol);
    std::vector<float> precomputedHist = SurfaceDetermination::computeHistogram(vol);

    auto getSlice = [&](int axis, int idx) {
        int w, h;
        if (axis == 0) { w = vol.dims[0]; h = vol.dims[1]; }
        else if (axis == 1) { w = vol.dims[0]; h = vol.dims[2]; }
        else               { w = vol.dims[1]; h = vol.dims[2]; }
        cv::Mat gray(h, w, CV_32F);
        for (int r = 0; r < h; ++r)
        for (int c = 0; c < w; ++c) {
            if (axis == 0)      gray.at<float>(r, c) = vol.at(c, r, idx);
            else if (axis == 1) gray.at<float>(r, c) = vol.at(c, idx, r);
            else                gray.at<float>(r, c) = vol.at(idx, c, r);
        }
        return gray;
    };
    cv::Mat sliceXY = getSlice(0, vol.dims[2]/2);

    // ── 3. 全参数交互控制 ────────────────────────────────────────
    std::string winName = "Parameter Controls";
    cv::namedWindow(winName, cv::WINDOW_AUTOSIZE);

    auto valToPos = [&](float val) { return static_cast<int>((val - vol.minVal) / (vol.maxVal - vol.minVal) * 1000); };
    auto posToVal = [&](int pos) { return vol.minVal + (pos / 1000.0f) * (vol.maxVal - vol.minVal); };

    int posISO = valToPos(params.isoValue);
    int posVoidThr = valToPos(params.background); 
    int posMinVol = 10; 
    int advDist = 20;

    cv::createTrackbar("Surface ISO", winName, &posISO, 1000);
    cv::createTrackbar("Void Threshold", winName, &posVoidThr, 1000);
    cv::createTrackbar("Min Vol(x10)", winName, &posMinVol, 1000);
    cv::createTrackbar("Adv: SearchDist", winName, &advDist, 100);

    AdvancedSurfaceParams advParams;
    VoidDetectionParams voidParams;

    while (true) {
        params.isoValue = posToVal(posISO);
        voidParams.grayMin = vol.minVal;
        voidParams.grayMax = posToVal(posVoidThr);
        voidParams.minVolumeMM3 = posMinVol / 10.0;
        advParams.enabled = true;
        advParams.normalSearchDistance = advDist / 10.0f;
        advParams.searchStep = 0.5f;

        params.background = voidParams.grayMax;
        cv::imshow(winName, SurfaceDetermination::visualizeHistogram(vol, precomputedHist, params));

        cv::Mat display, binary;
        cv::normalize(sliceXY, display, 0, 255, cv::NORM_MINMAX);
        display.convertTo(display, CV_8U);
        cv::cvtColor(display, display, cv::COLOR_GRAY2BGR);
        
        cv::threshold(sliceXY, binary, params.isoValue, 255, cv::THRESH_BINARY);
        binary.convertTo(binary, CV_8U);
        std::vector<std::vector<cv::Point>> contours;
        cv::findContours(binary, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
        SurfaceDetermination::refineContour2D(vol, 0, vol.dims[2]/2, contours, params, advParams);
        cv::drawContours(display, contours, -1, cv::Scalar(0, 255, 255), 1);

        for(int r=0; r<display.rows; r+=4) for(int c=0; c<display.cols; c+=4) {
            if(sliceXY.at<float>(r,c) <= voidParams.grayMax) {
                display.at<cv::Vec3b>(r,c) = cv::Vec3b(255, 0, 0);
            }
        }

        cv::imshow("XY Preview", display);
        int key = cv::waitKey(30);
        if (key == 13 || key == 27) break;
    }
    cv::destroyAllWindows();

    // ── 4. 批处理分析 ──────────────────────────────────────────
    printf("\nProcessing... Please wait.\n");
    std::vector<uint8_t> interiorMask = VoidAnalysis::createInteriorMaskVoxel(vol, params.isoValue);
    std::vector<uint8_t> candidateVoids = VoidAnalysis::extractCandidateVoids(vol, interiorMask, voidParams);
    std::vector<int> labelVolume;
    std::vector<VoidRegion> voids = VoidAnalysis::labelAndAnalyze(vol, candidateVoids, voidParams, labelVolume);

    // ── 4.1 导出结果 (Step 8 & 10) ──────────────────────────────
    VoidAnalysis::saveResults("output", voids, labelVolume, vol);

    /*
    
        // ── 4.2 生成 3D 效果 (Step 9.3) ─────────────────────────────
    printf("Generating 3D Void Mesh for visualization...\n");
    vtkSmartPointer<vtkPolyData> voidMesh = VoidAnalysis::generateVoidMesh(labelVolume, vol);
    
    auto colors = vtkSmartPointer<vtkNamedColors>::New();
    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(voidMesh);
    mapper->ScalarVisibilityOff();

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());

    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetWindowName("3D Void Visualization");

    auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);

    renderer->AddActor(actor);
    renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());

    printf("Starting 3D visualization. Use mouse to rotate. Close window to continue.\n");
    renderWindow->Render();
    interactor->Start();

    */
    // ── 5. 结果全轴预览 ────────────────────────────────────────
    std::string resWin = "Final Analysis Results";
    cv::namedWindow(resWin);
    int currentSlice = vol.dims[2] / 2;
    int viewAxis = 0;
    cv::createTrackbar("Axis", resWin, &viewAxis, 2);
    cv::createTrackbar("Slice", resWin, &currentSlice, vol.dims[2]-1);

    while (true) {
        int maxIdx = (viewAxis == 0) ? vol.dims[2]-1 : (viewAxis == 1 ? vol.dims[1]-1 : vol.dims[0]-1);
        cv::setTrackbarMax("Slice", resWin, maxIdx);
        if (currentSlice > maxIdx) currentSlice = maxIdx;

        int w, h;
        if (viewAxis == 0) { w = vol.dims[0]; h = vol.dims[1]; }
        else if (viewAxis == 1) { w = vol.dims[0]; h = vol.dims[2]; }
        else { w = vol.dims[1]; h = vol.dims[2]; }

        cv::Mat gray(h, w, CV_32F);
        for (int r = 0; r < h; ++r) for (int c = 0; c < w; ++c) {
            if (viewAxis == 0)      gray.at<float>(r, c) = vol.at(c, r, currentSlice);
            else if (viewAxis == 1) gray.at<float>(r, c) = vol.at(c, currentSlice, r);
            else                    gray.at<float>(r, c) = vol.at(currentSlice, c, r);
        }

        cv::Mat display;
        cv::normalize(gray, display, 0, 255, cv::NORM_MINMAX);
        display.convertTo(display, CV_8U);
        cv::cvtColor(display, display, cv::COLOR_GRAY2BGR);

        for (int r = 0; r < h; ++r) for (int c = 0; c < w; ++c) {
            size_t idx;
            if (viewAxis == 0)      idx = (size_t)c + (size_t)r * vol.dims[0] + (size_t)currentSlice * vol.dims[0] * vol.dims[1];
            else if (viewAxis == 1) idx = (size_t)c + (size_t)currentSlice * vol.dims[0] + (size_t)r * vol.dims[0] * vol.dims[1];
            else                    idx = (size_t)currentSlice + (size_t)c * vol.dims[0] + (size_t)r * vol.dims[0] * vol.dims[1];
            if (labelVolume[idx] > 0) {
                cv::Vec3b& p = display.at<cv::Vec3b>(r, c);
                p[2] = 255; p[0] /= 2; p[1] /= 2;
            }
        }
        cv::imshow(resWin, display);
        if (cv::waitKey(30) == 27) break;
    }
    
    
    

    return 0;
}
