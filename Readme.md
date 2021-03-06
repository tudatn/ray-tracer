## Introduction

A simple ray tracer implementation in C++ with openGL.

## Functions

The program implements Phong illumination model to calculate colors of intersection points. Additional options include:
- shadow
- reflection
- refraction
- diffusion
- super sampling 

Objects of the scene are spheres with different parameters (radius, position, ambient, specular, and diffuse). A chessboard which is modeled by a sphere is added to test the reflection.

## How to compile

```
make clean_object
make
./raycast [-d|-u] [0-5] [+s|+l|+c|+r|+f|+p]
```

## Sample rendered photos

![sample](https://user-images.githubusercontent.com/26643781/39025448-e3cae042-43fb-11e8-836a-91c51ca1aabe.png)
