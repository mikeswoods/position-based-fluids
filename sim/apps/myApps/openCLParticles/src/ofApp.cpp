/******
 
 This example updates 1M particles on the GPU using OpenCL
 The OpenCL kernel writes position data directly to a VBO stored in the OpenGL device memory
 so now data transfer between host and device during runtime
 
 
 Kernel based on Rui's ofxOpenCL particle example opencl particles 001b.zip
 at http://code.google.com/p/ruisource/
 *****/


#include "ofApp.h"
#include "MSAOpenCL.h"

#define NUM_PARTICLES (100000)


typedef struct{
    float4 pos;
    float4 vel;
    float mass;
    float dummy[3];		// need this to make sure the float2 vel is aligned to a 16 byte boundary
} Particle;


float4 mousePos;
float4 dimensions;


msa::OpenCL opencl;

// vector of Particles on host and corresponding clBuffer on device
msa::OpenCLBufferManagedT<Particle>	particles;

GLuint vbo;

//--------------------------------------------------------------
void ofApp::setup(){
    ofBackground(0, 0, 0);
    ofSetLogLevel(OF_LOG_VERBOSE);
    ofSetVerticalSync(false);
    
    opencl.setupFromOpenGL();
    
    // create vbo
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(Particle) * NUM_PARTICLES, 0, GL_DYNAMIC_COPY);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    // init host and CL buffers
    particles.initFromGLObject(vbo, NUM_PARTICLES);
    
    // init data
    
    for(int i=0; i<NUM_PARTICLES; i++) {
        Particle &p = particles[i];
        p.pos.x = ofRandomWidth();
        p.pos.y = ofRandomHeight();
        p.vel.x = 0;
        p.mass  = ofRandom(0.5, 1);
    }
    
    particles.writeToDevice();
    
    opencl.loadProgramFromFile("MSAOpenCL/Particle.cl");
    opencl.loadKernel("updateParticle");
    
    opencl.kernel("updateParticle")->setArg(0, particles);
    opencl.kernel("updateParticle")->setArg(1, mousePos);//.getPtr(), sizeof(float2));
    opencl.kernel("updateParticle")->setArg(2, dimensions);//.getPtr(), sizeof(float2));
    
    glPointSize(5);
}




//--------------------------------------------------------------
void ofApp::update(){
    
    mousePos.x = ofGetMouseX();
    mousePos.y = ofGetMouseY();
    dimensions.x = ofGetWidth();
    dimensions.y = ofGetHeight();
    
    opencl.kernel("updateParticle")->setArg(1, mousePos);//.getPtr(), sizeof(float2));
    opencl.kernel("updateParticle")->setArg(2, dimensions);//.getPtr(), sizeof(float2) );
    glFlush();
    
    opencl.kernel("updateParticle")->run1D(NUM_PARTICLES);
}


//--------------------------------------------------------------
void ofApp::draw(){
    opencl.finish();
    
    glColor3f(1.0f, 1.0f, 1.0f);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    
    opencl.finish();
    
    //glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_FLOAT, 0, 0);
    glDrawArrays(GL_POINTS, 0, NUM_PARTICLES);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    
    glColor3f(1, 1, 1);
    string info = "fps: " + ofToString(ofGetFrameRate()) + "\nnumber of particles: " + ofToString(NUM_PARTICLES);
    ofDrawBitmapString(info, 20, 20);
}
