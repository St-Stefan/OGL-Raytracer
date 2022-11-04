#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#include "texture.h"
#include <math.h>
#include <numbers>
#ifdef NDEBUG
#include <omp.h>
#endif

int depth = 5;
const int numRays = 8; 
//Implementing the recursive ray-tracer 
//Source: Chapter 4.8, Fundamentals of Computer Graphics, Fourth Edition.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {
        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        drawRay(ray, Lo);

         Ray reflection = computeReflectionRay(ray, hitInfo);


        if (features.enableRecursive) {
          if (rayDepth >= depth || hitInfo.material.ks == glm::vec3 { 0.0f }) {

          } 
          else {
              if (features.extra.enableGlossyReflection) {
                  glm::vec3 color { 0.0f };
                  auto W = glm::abs(hitInfo.normal.x) > 0.99 ? glm::vec3 { 0, 1, 0 } : glm::vec3 { 1, 0, 0 };
                  glm::vec3 u_vector = glm::normalize(glm::cross(hitInfo.normal, W));
                  glm::vec3 v_vector = glm::normalize(glm::cross(reflection.direction, u_vector));

                  for (int i = 0; i < numRays; i++) {
                      for (int j = 0; j < numRays; j++) {
                          float u = -1 / (2.0f * hitInfo.material.shininess) + (float)std::rand() / RAND_MAX * 1 / hitInfo.material.shininess;
                          float v = -1 / (2.0f * hitInfo.material.shininess) + (float)std::rand() / RAND_MAX * 1 / hitInfo.material.shininess;

                          glm::vec3 r_dash = reflection.direction + u * u_vector + v * v_vector;

                          float weight = glm::pow(glm::dot(reflection.direction, glm::normalize(r_dash)), hitInfo.material.shininess);

                          color += getFinalColor(scene, bvh, Ray { reflection.origin, r_dash }, features, rayDepth + 1) * weight;
                      }
                  }
                  color /= (numRays * numRays);
                  return Lo += hitInfo.material.ks * color;
              }

              Lo += hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
          }
        }


        if (features.extra.enableTransparency) {
            float offset = 0.001f; // offset to remove acne from transparent surface
            Ray transparency;
            transparency.direction = ray.direction; // the new ray goes through the surface, in the same direction
            transparency.origin = ray.origin + ray.direction * ray.t + offset * ray.direction; // offset is added to the origin of the new ray
            if (rayDepth < 10 && hitInfo.material.transparency != 1) { // check to make sure ray does not pass through opaque materials
                Lo = hitInfo.material.transparency * Lo + (1 - hitInfo.material.transparency) * getFinalColor(scene, bvh, transparency, features, rayDepth + 1);
            }
        }
     
        // Visual Debug: Draw a ray with a color which is the returned value from computeLightContribution
  
        // Set the color of the pixel to white if the ray hits.
        return Lo;

    } else {
        // Draw a red debug ray if the ray missed.
        // Set the color of the pixel to black if the ray misses.
        
         drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
         return glm::vec3(0.0f);
        
    }
}

glm::vec3 Aux(const Trackball& camera, Ray cameraRay, const BvhInterface bvh, const Features& features, const Scene& scene, Screen& screen)
{
    float focusDistance = features.extra.FocusDistance;
    float aperture = features.extra.Aperture;
    int samples = features.extra.Samples;
    glm::vec3 pixelColor = { 0, 0, 0 };
    float g1 = ((float)rand() / (float)(RAND_MAX / 2)) - 1;
    float g2 = ((float)rand() / (float)(RAND_MAX / 2)) - 1;
    float g3 = ((float)rand() / (float)(RAND_MAX / 2)) - 1;
    glm::vec3 perpendic = glm::cross(camera.forward(), glm::vec3 { g1, g2, g3 });
    perpendic = glm::normalize(perpendic) * aperture;
    glm::vec3 focusPoint = cameraRay.origin + cameraRay.direction * focusDistance;
    for (int i = 0; i < samples; i++) {
        Ray G;
        G.origin = cameraRay.origin + perpendic;
        G.direction = glm::normalize(focusPoint - G.origin);
        pixelColor += getFinalColor(scene, bvh, G, features);
    }
    pixelColor = pixelColor / float(samples);
    return pixelColor;
}

std::vector<float> kernelCreation(int filterSize)
{
    std::vector<float> kernel((2 * filterSize + 1) * (2 * filterSize + 1));
    // initialising standard deviation to 1.0
    float sigma = 1.0f;
    float r, s = 2.0f * sigma * sigma;

    // sum is for normalization
    float sum = 0.0;

    // generating 5x5 kernel
    for (int x = -filterSize; x <= filterSize; x++) {
        for (int y = -filterSize; y <= filterSize; y++) {
            int x2 = x + filterSize;
            int y2 = y + filterSize;
            r = (float)sqrt(x * x + y * y);
            kernel[(2 * filterSize + 1 - 1 - y2) * (2 * filterSize + 1) + x2] = (exp(-(r * r) / s)) / (glm::pi<float>() * s);
            sum += kernel[(2 * filterSize + 1 - 1 - y2) * (2 * filterSize + 1) + x2];
        }
    }

    // normalising the Kernel
    for (int x = -filterSize; x <= filterSize; ++x) {
        for (int y = -filterSize; y <= filterSize; ++y) {
            int x2 = x + filterSize;
            int y2 = y + filterSize;
            kernel[(2 * filterSize + 1 - 1 - y2) * (2 * filterSize + 1) + x2] /= sum;
        }
    }
    return kernel;
}

glm::vec3 pixelAt(int i, int j, Screen& screen)
{
    if (i < 0 || i >= screen.resolution().x || j < 0 || j >= screen.resolution().y) {
        return glm::vec3 { 0.0f };
    }
    return screen.pixels()[screen.indexAt(i, j)];
}

glm::vec3 boxFilter(int i, int j, int filterSize, Screen& screen, float scale, std::vector<float> kernel)
{
    filterSize = std::max(1, filterSize);
    glm::vec3 result { 0.0f };
    for (int x = -filterSize; x < filterSize + 1; ++x) {
        for (int y = -filterSize; y < filterSize + 1; ++y) {
            int x2 = x + filterSize;
            int y2 = y + filterSize;
            result += kernel[(2 * filterSize + 1 - 1 - y2) * (2 * filterSize + 1) + x2] * pixelAt(i + x, j + y, screen);
        }
    }
    result *= scale;
    return result;
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    float focusDistance = features.extra.FocusDistance;
    float aperture = features.extra.Aperture;
    int samples = features.extra.Samples;

    int width = windowResolution.x;
    int height = windowResolution.y;
    glm::vec3 inital { 0.0f };
    std::vector<glm::vec3> pixels(width * height);

    Screen mask = Screen(screen.resolution());

    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            glm::vec3 pixelColor = { 0, 0, 0 };
            if (features.extra.enableDepthOfField) {
                pixelColor = Aux(camera, cameraRay, bvh, features, scene, screen);
                screen.setPixel(x, y, pixelColor);
            } else if (features.extra.enableMultipleRaysPerPixel) {
                glm::vec3 colorNotDivided = glm::vec3 { 0.f };
                colorNotDivided += getFinalColor(scene, bvh, cameraRay, features);
                if (features.extra.enableMultipleRaysPerPixel) {
                    int samplesHere = features.extra.samplesPerPixel;
                    for (int i = 1; i < samplesHere; ++i) {
                        float angleInRadians = 2.f * std::numbers::pi_v<float> * rand() / (float)RAND_MAX;
                        float r = 1.f / float(windowResolution.x * 2) * 2.f * rand() / (float)RAND_MAX;
                        glm::vec2 offset = glm::vec2 { std::cos(angleInRadians) * r, std::sin(angleInRadians) * r };
                        glm::vec2 localPixelPos = normalizedPixelPos + offset;
                        const Ray cameraRay1 = camera.generateRay(localPixelPos);
                        colorNotDivided += getFinalColor(scene, bvh, cameraRay1, features);
                    }

                    screen.setPixel(x, y, colorNotDivided / float(samplesHere));
                } else {
                    screen.setPixel(x, y, colorNotDivided);
                }
            } else {
                pixelColor = getFinalColor(scene, bvh, cameraRay, features);
                    screen.setPixel(x, y, pixelColor);
            }
        }
    }

    if (features.extra.enableBloomEffect) {
        float scale = 2.0f;
        int filterSize = 8;
        float threshold = 0.6f;
        std::vector<float> kernel = kernelCreation(filterSize);
        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                glm::vec3 color = screen.pixels()[screen.indexAt(x, y)];
                if (color.x > threshold && color.y > threshold && color.z > threshold) {
                    pixels[screen.indexAt(x, y)] = color;
                } else {
                    pixels[screen.indexAt(x, y)] = glm::vec3 { 0.0f };
                }
            }
        }

        mask.pixels() = pixels;

        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                pixels[screen.indexAt(x, y)] = boxFilter(x, y, filterSize, mask, scale, kernel);
            }
        }

        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                if (pixels[screen.indexAt(x, y)] != glm::vec3 { 0.0f }) {
                    screen.setPixel(x, y, (screen.pixels()[screen.indexAt(x, y)] + pixels[screen.indexAt(x, y)]));
                }
            }
        }
    }
}