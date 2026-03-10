#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

int main() {
    int dimX = 64, dimY = 64, dimZ = 64;
    std::vector<float> data(dimX * dimY * dimZ, 0.0f);
    
    float cx = dimX / 2.0f;
    float cy = dimY / 2.0f;
    float cz = dimZ / 2.0f;
    float radius = 20.0f;
    
    for (int z = 0; z < dimZ; ++z) {
        for (int y = 0; y < dimY; ++y) {
            for (int x = 0; x < dimX; ++x) {
                float dx = x - cx;
                float dy = y - cy;
                float dz = z - cz;
                float dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                // Background is ~0, Sphere is ~255
                float val = 0.0f;
                if (dist < radius) {
                    val = 255.0f;
                } else if (dist < radius + 2.0f) {
                    // Smooth transition
                    val = 255.0f * (1.0f - (dist - radius) / 2.0f);
                }
                
                data[x + y * dimX + z * dimX * dimY] = val;
            }
        }
    }
    
    // Write RAW file
    std::ofstream rawFile("test_volume.raw", std::ios::binary);
    rawFile.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
    rawFile.close();
    
    // Write MHD file
    std::ofstream mhdFile("test_volume.mhd");
    mhdFile << "ObjectType = Image\n";
    mhdFile << "NDims = 3\n";
    mhdFile << "BinaryData = True\n";
    mhdFile << "BinaryDataByteOrderMSB = False\n";
    mhdFile << "CompressedData = False\n";
    mhdFile << "TransformMatrix = 1 0 0 0 1 0 0 0 1\n";
    mhdFile << "Offset = 0 0 0\n";
    mhdFile << "CenterOfRotation = 0 0 0\n";
    mhdFile << "ElementSpacing = 1 1 1\n";
    mhdFile << "DimSize = " << dimX << " " << dimY << " " << dimZ << "\n";
    mhdFile << "ElementType = MET_FLOAT\n";
    mhdFile << "ElementDataFile = test_volume.raw\n";
    mhdFile.close();
    
    std::cout << "Generated test_volume.mhd and test_volume.raw" << std::endl;
    return 0;
}
