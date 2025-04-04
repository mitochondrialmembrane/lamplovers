#version 330 core
out vec4 fragColor;

// Additional information for lighting
in vec4 normal_worldSpace;
in vec4 position_worldSpace;

uniform int wire = 0;
uniform float red = 0;
uniform float green = 0;
uniform float blue = 0;
uniform float alpha = 1;
uniform int u_renderMode;

void main() {
    if (u_renderMode == 1) { // Point Cloud Shading
        fragColor = vec4(red, green, blue, alpha);
        // Simple color-based point rendering
        //float dist = length(gl_PointCoord - vec2(0.5)); // Circular effect
        //if (dist > 0.5) discard; // Make points round

        //fragColor = vec4(red, green, blue, 1.0);
    }
    else {
        if (wire == 1) {
            fragColor = vec4(0.0, 0.0, 0.0, 1);
            return;
        }
        vec4 lightPos   = vec4(-2.0, 2.0, -3.0 , 1.0);
        vec3 lightColor = vec3(1.0f, alpha, 0.0f);
        vec4 lightDir   = normalize(-lightPos + position_worldSpace);
        float c = clamp(dot(-normal_worldSpace, lightDir), 0, 1);

        fragColor = vec4(red * c * lightColor[0], green * c * lightColor[0], blue * c * lightColor[0], alpha);
    }
}
