#include "texture.h"
#include <framework/image.h>
#include <glm/geometric.hpp>
#include <shading.h>


//Source: Chapter 11.1, Fundamentals of Computer Graphics
glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
  
    int i = std::roundf(texCoord.x * image.width - 0.5f);
    int j = std::roundf((1-texCoord.y) * image.height - 0.5f);
    int index = j * image.width + i;
    return image.pixels[index];
 
}

/* std::pair<glm::vec2, int> cubeMapLookUp(glm::vec3 direction)
{
    glm::vec3 vAbs = glm::abs(direction);
    float faceIndex = 0.0f;
    float ma;
    glm::vec2 uv;
    if (vAbs.z >= vAbs.x && vAbs.z >= vAbs.y) {
        faceIndex = direction.z < 0.0 ? 5.0 : 4.0;
        ma = 0.5 / vAbs.z;
        uv = glm::vec2(direction.z < 0.0 ? -direction.x : direction.x, -direction.y);
    } else if (vAbs.y >= vAbs.x) {
        faceIndex = direction.y < 0.0 ? 3.0 : 2.0;
        ma = 0.5 / vAbs.y;
        uv = glm::vec2(direction.x, direction.y < 0.0 ? -direction.z : direction.z);
    } else {
        faceIndex = direction.x < 0.0 ? 1.0 : 0.0;
        ma = 0.5 / vAbs.x;
        uv = glm::vec2(direction.x < 0.0 ? direction.z : -direction.z, -direction.y);
    }
    return { uv * ma + 0.5f, faceIndex };
}*/


std::pair<std::pair<float,float>, int> cubeMapLookUp(float x, float y, float z)
{
    float absX = fabs(x);
    float absY = fabs(y);
    float absZ = fabs(z);
    int index = 0; 
    float u = 0.0f; 
    float v = 0.0f;

    int isXPositive = x > 0 ? 1 : 0;
    int isYPositive = y > 0 ? 1 : 0;
    int isZPositive = z > 0 ? 1 : 0;

    float maxAxis, uc, vc;

    // POSITIVE X
    if (isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from +z to -z
        // v (0 to 1) goes from -y to +y
        maxAxis = absX;
        uc = -z;
        vc = y;
        index = 0;
    }
    // NEGATIVE X
    if (!isXPositive && absX >= absY && absX >= absZ) {
        // u (0 to 1) goes from -z to +z
        // v (0 to 1) goes from -y to +y
        maxAxis = absX;
        uc = z;
        vc = y;
        index = 1;
    }
    // POSITIVE Y
    if (isYPositive && absY >= absX && absY >= absZ) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from +z to -z
        maxAxis = absY;
        uc = x;
        vc = -z;
        index = 2;
    }
    // NEGATIVE Y
    if (!isYPositive && absY >= absX && absY >= absZ) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from -z to +z
        maxAxis = absY;
        uc = x;
        vc = z;
        index = 3;
    }
    // POSITIVE Z
    if (isZPositive && absZ >= absX && absZ >= absY) {
        // u (0 to 1) goes from -x to +x
        // v (0 to 1) goes from -y to +y
        maxAxis = absZ;
        uc = x;
        vc = y;
        index = 4;
    }
    // NEGATIVE Z
    if (!isZPositive && absZ >= absX && absZ >= absY) {
        // u (0 to 1) goes from +x to -x
        // v (0 to 1) goes from -y to +y
        maxAxis = absZ;
        uc = -x;
        vc = y;
        index = 5;
    }

    // Convert range from -1 to 1 to 0 to 1
    u = 0.5f * (uc / maxAxis + 1.0f);
    v = 0.5f * (vc / maxAxis + 1.0f);

    return {{ u, v }, index};
}

glm::vec3 environmentMapping(const std::vector<Image>& images,const HitInfo& hitInfo,
        const Ray& ray, const Features& features)
{  
    std::pair<std::pair<float, float>, int> result = cubeMapLookUp(ray.direction.x, ray.direction.y, ray.direction.z);
    return acquireTexel(images[result.second], glm::vec2{result.first.first, result.first.second} , features);
}