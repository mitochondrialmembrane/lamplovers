#version 330 core

layout(location = 0) in vec3 position; // Position of the vertex
layout(location = 1) in float temp;
// layout(location = 1) in vec3 normal;   // Normal of the vertex

uniform mat4 proj;
uniform mat4 view;
uniform mat4 model;

uniform mat3 inverseTransposeModel;
uniform float pointSize = 1.0;
uniform int u_renderMode;

// out vec4 normal_worldSpace;
out vec4 position_worldSpace;
out float vTemp;

void main() {
    if (u_renderMode == 1) {
        gl_PointSize = pointSize;
    }
    // else {
    //     normal_worldSpace = vec4(normalize(inverseTransposeModel * normal), 0);
    // }
    position_worldSpace = vec4(position, 1.0);

    gl_Position = proj * view * model * vec4(position, 1.0);

    vTemp = temp;
}
