#version 150

uniform float particleRadius;
uniform vec3  cameraPosition;
uniform mat4  modelViewProjectionMatrix;

out vec4 color;
in vec4 position;
in vec4 normal;

void main()
{
    float minPointScale = 1.0e-6;
    float maxPointScale = particleRadius * 2.0;
    float maxDistance   = 100.0;
    
    gl_Position = modelViewProjectionMatrix * position;
    
    float cameraDist = distance(position.xyz, cameraPosition);
    float pointScale = clamp(1.0 - (cameraDist / maxDistance), minPointScale, maxPointScale);
    
    float wScale = position.w;
    if (wScale == 0.0) {
        wScale = 0.0001;
    }
    
    gl_PointSize = (wScale * pointScale) * particleRadius;
    
    // Waterish color:
    color = vec4(0.251, 0.643, 0.875, 0.8);
}
