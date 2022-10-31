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
    int j = std::roundf(texCoord.y * image.height - 0.5f);
    int index = j * image.width + i;
    return image.pixels[index];
 
}

std::pair<glm::vec2, int> cubeMapLookUp(glm::vec3 direction)
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
}

glm::vec3 environmentMapping(const std::vector<Image>& images,const HitInfo& hitInfo,
        const Ray& ray, const Features& features, bool hit)
{
 
    Ray reflection = computeReflectionRay(ray, hitInfo);
    std::pair<glm::vec2, int> result;

    if (hit == false) {
         result = cubeMapLookUp(ray.direction);
    } else {
        result = cubeMapLookUp(reflection.direction);
    }
    
    return acquireTexel(images[result.second], result.first, features);
}