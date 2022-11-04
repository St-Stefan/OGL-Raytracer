#include "texture.h"
#include <framework/image.h>


//Source: Chapter 11.1, Fundamentals of Computer Graphics
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

//Source: Chapter 11.3.2, Fundamentals of Computer Graphics 
glm::vec3 bilinearInterpolation(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    float u_p = texCoord.x * image.width - 0.5;
    float v_p = (1 - texCoord.y) * image.height - 0.5;
    float iu0 = std::floorf(u_p);
    float iu1 = iu0 + 1;
    float iv0 = std::floorf(v_p);
    float iv1 = iv0 + 1;
    float a_u = (iu1 - u_p);
    float b_u = 1 - a_u;
    float a_v = (iv1 - v_p);
    float b_v = 1 - a_v;
    return a_u * a_v * acquireTexelClampMode(iu0, iv0, image)
        + a_u * b_v *  acquireTexelClampMode(iu0, iv1, image)
        + b_u * a_v *  acquireTexelClampMode(iu1, iv0, image)
        + b_u * b_v *    acquireTexelClampMode(iu1, iv1, image);
}