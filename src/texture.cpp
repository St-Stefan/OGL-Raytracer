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
    int j = std::roundf(texCoord.y * image.height - 0.5f);
    int index = j * image.width + i;
    return image.pixels[index];
    if (index >= image.pixels.size()) {
        index = image.pixels.size() - 1;
    } 
    if (index < 0) {
        index = 0;
    }
    return image.pixels[index];
}