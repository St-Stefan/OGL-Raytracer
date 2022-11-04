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

std::pair<std::pair<float,float>, int> cubeMapLookUp(float x_direction, float y_direction, float z_direction)
{
    float absolute_x = std::abs(x_direction);
    float absolute_y = std::abs(y_direction);
    float absolute_z = std::abs(z_direction);
    int index = 0; 
    float u = 0.0f; 
    float v = 0.0f;

    int isXPositive = x_direction > 0 ? 1 : 0;
    int isYPositive = y_direction > 0 ? 1 : 0;
    int isZPositive = z_direction > 0 ? 1 : 0;

    float axis = 0.0f;
    float u_coordinate = 0.0f;
    float v_coordinate = 0.0f;

    // POSITIVE X
    if (isXPositive && absolute_x >= absolute_y && absolute_x >= absolute_z) {
       
        axis = absolute_x;
        u = -z_direction;
        v = y_direction;
        index = 0;
    }
    // NEGATIVE X
    if (!isXPositive && absolute_x >= absolute_y && absolute_x >= absolute_z) {
        
        axis = absolute_x;
        u = z_direction;
        v = y_direction;
        index = 1;
    }
    // POSITIVE Y
    if (isYPositive && absolute_y >= absolute_x && absolute_y >= absolute_z) {
        
        axis = absolute_y;
        u = x_direction;
        v = -z_direction;
        index = 2;
    }
    // NEGATIVE Y
    if (!isYPositive && absolute_y >= absolute_x && absolute_y >= absolute_z) {
        
        axis = absolute_y;
        u = x_direction;
        v = z_direction;
        index = 3;
    }
    // POSITIVE Z
    if (isZPositive && absolute_z >= absolute_x && absolute_z >= absolute_y) {
        
        axis = absolute_z;
        u = x_direction;
        v = y_direction;
        index = 4;
    }
    // NEGATIVE Z
    if (!isZPositive && absolute_z >= absolute_x && absolute_z >= absolute_y) {
       
        axis = absolute_z;
        u = -x_direction;
        v = y_direction;
        index = 5;
    }

    // Convert range from -1 to 1 to 0 to 1
    u_coordinate = 0.5f * (u / axis + 1.0f);
    v_coordinate = 0.5f * (v / axis + 1.0f);

    return {{ u_coordinate, v_coordinate }, index};
}

glm::vec3 environmentMapping(const std::vector<Image>& images,
        const Ray& ray, const Features& features)
{  
    std::pair<std::pair<float, float>, int> result = cubeMapLookUp(ray.direction.x, ray.direction.y, ray.direction.z);
    return acquireTexel(images[result.second], glm::vec2{result.first.first, result.first.second} , features);
}