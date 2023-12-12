#include <emscripten.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

typedef unsigned char byte;
typedef unsigned int uint32;
typedef unsigned long long uint64;
typedef float float32;

typedef struct {
    unsigned char r, g, b, a;
} Pixel;

typedef enum {
    PROJECTION_DIMETRIC,
    PROJECTION_TRIMETRIC
} Projection;

typedef enum {
    DIRECTION_CW,
    DIRECTION_CCW
} Direction;

bool isTransparent(Pixel p) {
    return p.a == 0;
}

Pixel interpolatePixel(Pixel left, Pixel right) {
    Pixel result;

    result.r = (left.r + right.r) / 2;
    result.g = (left.g + right.g) / 2;
    result.b = (left.b + right.b) / 2;
    result.a = (left.a + right.a) / 2;

    return result;
}

float getBrightness(Pixel p) {
    return (0.299*p.r + 0.587*p.g + 0.114*p.b) / 255;
}

float getIntensity(Pixel p) {
    return (p.r + p.g + p.b) / (3.0 * 255);
}

float getSaturation(Pixel p) {
    float r = p.r / 255.0;
    float g = p.g / 255.0;
    float b = p.b / 255.0;

    float min = fmin(fmin(r, g), b);
    float max = fmax(fmax(r, g), b);
    float delta = max - min;

    if(max == 0) {
        return 0;
    }

    return delta / max;
}

// Переводим RGB в HSV и возвращаем высоту в соответствии с направлением обхода и начальным цветом
float RGBtoHSV_direction(Pixel p, float startHue, bool clockwise) {
    float r = p.r / 255.0f;
    float g = p.g / 255.0f;
    float b = p.b / 255.0f;

    float max = fmaxf(fmaxf(r, g), b);
    float min = fminf(fminf(r, g), b);
    float delta = max - min;

    float hue;
    if (delta < 0.00001f)
        hue = 0;
    else if (max == r)
        hue = (g - b) / delta;
    else if (max == g)
        hue = 2 + (b - r) / delta;
    else
        hue = 4 + (r - g) / delta;

    hue *= 60;
    if (hue < 0)
        hue += 360;

    // Циклически сдвигаем, чтобы начальный цвет был на нуле
    hue = fmodf(hue - startHue + 360, 360);

    // Возвращаем инвертированную высоту, если направление по часовой стрелке
    return clockwise ? (360 - hue) / 360 : hue / 360;
}

// Критерии высоты
float getRedCW(Pixel p) {
    return RGBtoHSV_direction(p, 0, true);
}

float getRedCCW(Pixel p) {
    return RGBtoHSV_direction(p, 0, false);
}

float getGreenCW(Pixel p) {
    return RGBtoHSV_direction(p, 120, true);
}

float getGreenCCW(Pixel p) {
    return RGBtoHSV_direction(p, 120, false);
}

float getBlueCW(Pixel p) {
    return RGBtoHSV_direction(p, 240, true);
}

float getBlueCCW(Pixel p) {
    return RGBtoHSV_direction(p, 240, false);
}


void apply3D(Pixel* image, int width, int height, float amount, const char* criterion,Projection projection, Direction direction) {
	
	// Устанавливаем все пиксели полностью непрозрачными
    for(int y = 0; y < height; ++y) {
        for(int x = 0; x < width; ++x) {
            image[y * width + x].a = 255;
        }
    }
	
    // Копируем исходное изображение во временное хранилище
    Pixel* tempImage = malloc(width * height * sizeof(Pixel));
    memcpy(tempImage, image, width * height * sizeof(Pixel));

    // Устанавливаем исходное изображение в полностью прозрачное
    memset(image, 0, width * height * sizeof(Pixel));

    float diagonal = sqrt(width * width + height * height);
    int corners[4][2] = {{0, 0}, {width - 1, 0}, {0, height - 1}, {width - 1, height - 1}};
    int minX = width, maxX = 0, maxY = 0;

    // Находим максимальные и минимальные координаты после проекции углов изображения
    for(int i = 0; i < 4; i++) {
        int x = corners[i][0];
        int y = corners[i][1];
        int newX, newY;

        switch(projection) {
            case PROJECTION_DIMETRIC:
                newX = direction == DIRECTION_CW ? (int)(width/2.0 + ((width * (x - y))/(2.0 * diagonal))) : 
                                                    (int)(width/2.0 + ((width * (y - x))/(2.0 * diagonal)));
                newY = (int)(height/2.0 + ((height * (x + y)) / (4.0 * diagonal)));
                break;
            case PROJECTION_TRIMETRIC:
                newX = direction == DIRECTION_CW ? (int)(width/2.0 + ((width * (2*x - y)) / (2.0 * diagonal))) : 
                                                    (int)(width/2.0 + ((width * (2*y - x)) / (2.0 * diagonal)));
                newY = (int)(height/2.0 + ((height * (x + 2*y)) / (4.0 * diagonal)));
                break;
        }

        minX = newX < minX ? newX : minX;
        maxX = newX > maxX ? newX : maxX;
        maxY = newY > maxY ? newY : maxY;
    }

    // Применяем проекцию с учетом найденных максимальных и минимальных значений
    for(int y = 0; y < height; ++y) {
        for(int x = 0; x < width; ++x) {
            int newX, newY;

            switch(projection) {
                case PROJECTION_DIMETRIC:
                    newX = direction == DIRECTION_CW ? (int)(width/2.0 + ((width * (x - y))/(2.0 * diagonal))) : 
                                                    (int)(width/2.0 + ((width * (y - x))/(2.0 * diagonal)));
                    newY = (int)(height/2.0 + ((height * (x + y)) / (4.0 * diagonal)));
                    break;
                case PROJECTION_TRIMETRIC:
                    newX = direction == DIRECTION_CW ? (int)(width/2.0 + ((width * (2*x - y)) / (2.0 * diagonal))) : 
                                                    (int)(width/2.0 + ((width * (2*y - x)) / (2.0 * diagonal)));
                    newY = (int)(height/2.0 + ((height * (x + 2*y)) / (4.0 * diagonal)));
                    break;
            }

            // Нормализация координат по x и корректировка y
            newX = ((newX - minX) * width) / (maxX - minX);
            newY = height - (maxY - newY); // Смещение изображения вниз

            // Проверка границ
            if(newX >= 0 && newX < width && newY >= 0 && newY < height) {
                // Копирование пикселя в новую позицию
                image[newY * width + newX] = tempImage[y * width + x];
            }
        }
    }

    // Интерполяция прозрачных пикселей
    for(int y = 0; y < height; ++y) {
        for(int x = 0; x < width; ++x) {
            if(isTransparent(image[y * width + x])) {
                int leftX = x - 1;
                while(leftX >= 0 && isTransparent(image[y * width + leftX])) {
                    leftX--;
                }

                int rightX = x + 1;
                while(rightX < width && isTransparent(image[y * width + rightX])) {
                    rightX++;
                }

                if(leftX >= 0 && rightX < width) {
                    image[y * width + x] = interpolatePixel(image[y * width + leftX], image[y * width + rightX]);
                }
            }
        }
    }
	
	// Освобождаем память временного изображения
    free(tempImage);
	
	////////////////Гистограммы//////////////////////
	Pixel* histogramTempImage = malloc(width * height * sizeof(Pixel));
	memcpy(histogramTempImage, image, width * height * sizeof(Pixel));
	
	float maxCriterion = 0;
    for(int y = 0; y < height; ++y) {
        for(int x = 0; x < width; ++x) {
            if(!isTransparent(image[y * width + x])) {
                float currentCriterion = 0;
                if(strcmp(criterion, "brightness") == 0) {
					currentCriterion = getBrightness(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "intensity") == 0) {
					currentCriterion = getIntensity(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "saturation") == 0) {
					currentCriterion = getSaturation(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "brightness_invert") == 0) {
					currentCriterion = 1 - getBrightness(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "intensity_invert") == 0) {
					currentCriterion = 1 - getIntensity(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "saturation_invert") == 0) {
					currentCriterion = 1 - getSaturation(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "RedCW") == 0) {
					currentCriterion = getRedCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "RedCCW") == 0) {
					currentCriterion = getRedCCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "GreenCW") == 0) {
					currentCriterion = getGreenCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "GreenCCW") == 0) {
					currentCriterion = getGreenCCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "BlueCW") == 0) {
					currentCriterion = getBlueCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "BlueCCW") == 0) {
					currentCriterion = getBlueCCW(histogramTempImage[y * width + x]);
				}

                maxCriterion = fmax(maxCriterion, currentCriterion);
            }
        }
    }

	for(int y = 0; y < height; ++y) {
		for(int x = 0; x < width; ++x) {
			if(!isTransparent(histogramTempImage[y * width + x])) {
				float currentCriterion = 0;
				if(strcmp(criterion, "brightness") == 0) {
					currentCriterion = getBrightness(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "intensity") == 0) {
					currentCriterion = getIntensity(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "saturation") == 0) {
					currentCriterion = getSaturation(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "brightness_invert") == 0) {
					currentCriterion = 1 - getBrightness(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "intensity_invert") == 0) {
					currentCriterion = 1 - getIntensity(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "saturation_invert") == 0) {
					currentCriterion = 1 - getSaturation(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "RedCW") == 0) {
					currentCriterion = getRedCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "RedCCW") == 0) {
					currentCriterion = getRedCCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "GreenCW") == 0) {
					currentCriterion = getGreenCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "GreenCCW") == 0) {
					currentCriterion = getGreenCCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "BlueCW") == 0) {
					currentCriterion = getBlueCW(histogramTempImage[y * width + x]);
				} else if(strcmp(criterion, "BlueCCW") == 0) {
					currentCriterion = getBlueCCW(histogramTempImage[y * width + x]);
				}

				int lineHeight = (int)((currentCriterion * height * amount) / (maxCriterion * 100));
				int startY = fmax(0, y - lineHeight);
				for(int lineY = startY; lineY < y; ++lineY) {
					float gradient = (float)(y - lineY) / lineHeight;
					image[lineY * width + x].r = histogramTempImage[y * width + x].r * gradient;
					image[lineY * width + x].g = histogramTempImage[y * width + x].g * gradient;
					image[lineY * width + x].b = histogramTempImage[y * width + x].b * gradient;
					image[lineY * width + x].a = histogramTempImage[y * width + x].a;
				}
			}
		}
	}

	free(histogramTempImage);
    
}

EMSCRIPTEN_KEEPALIVE
byte* wasmAlloc(uint32 width, uint32 height) {
  return malloc(width * height * 4);
}

EMSCRIPTEN_KEEPALIVE
void wasmFree(byte* p) {
  free(p);
}

EMSCRIPTEN_KEEPALIVE
void wasmProcess(Pixel* inputData, uint32 size, int width, int height, float amount, const char* criterion, const char* projection, const char* direction) {
	
	if (strcmp(projection, "dimetric") == 0) {
		if (strcmp(direction, "cw") == 0) {
			apply3D(inputData, width, height, amount, criterion, PROJECTION_DIMETRIC, DIRECTION_CW);
		} else  {
			apply3D(inputData, width, height, amount, criterion, PROJECTION_DIMETRIC, DIRECTION_CCW);
		} 
	} else if (strcmp(projection, "trimetric") == 0) {
		if (strcmp(direction, "cw") == 0) {
			apply3D(inputData, width, height, amount, criterion, PROJECTION_TRIMETRIC, DIRECTION_CW);
		} else  {
			apply3D(inputData, width, height, amount, criterion, PROJECTION_TRIMETRIC, DIRECTION_CCW);
		}
	} 
	
}



