#include "texture.h"
#include <framework/image.h>
#include <glm/geometric.hpp>
#include <shading.h>


// Source: Chapter 11.1, Fundamentals of Computer Graphics
glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    // TODO: implement this function.
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    int i = std::roundf(texCoord.x * image.width - 0.5f);
    int j = std::roundf((1 - texCoord.y) * image.height - 0.5f);
    int index = j * image.width + i;
    return acquireTexelClampMode(i, j, image);
}

glm::vec3 acquireTexelClampMode(int i, int j, const Image& image)
{
    if (j * image.width + i >= image.pixels.size()) {
        return image.pixels[image.pixels.size() - 1];
    }
    if (j * image.width + i < 0) {
        return image.pixels[0];
    }
    return image.pixels[j * image.width + i];
}

// Source: Chapter 11.3.2, Fundamentals of Computer Graphics
glm::vec3 bilinearInterpolation(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    float uinnerPoint = texCoord.x * image.width - 0.5;
    float vinnerPoint = (1 - texCoord.y) * image.height - 0.5;
    float u0 = std::floorf(uinnerPoint);
    float u1 = u0 + 1;
    float v0 = std::floorf(vinnerPoint);
    float v1 = v0 + 1;
    float alphaU = (u1 - uinnerPoint);
    float betaU = 1 - alphaU;
    float alphaV = (v1 - vinnerPoint);
    float betaV = 1 - alphaV;
    return alphaU * alphaV * acquireTexelClampMode(u0, v0, image)
        + alphaU * betaV * acquireTexelClampMode(u0, v1, image)
        + betaU * alphaV * acquireTexelClampMode(u1, v0, image)
        + betaU * betaV * acquireTexelClampMode(u1, v1, image);
}