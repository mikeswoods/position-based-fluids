/*******************************************************************************
 * ofApp
 * - The OpenFrameworks application entry point
 *
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#include <memory>
#include <iostream>
#include "ofApp.h"
#include "Definitions.h"

/******************************************************************************/

using namespace std;

/******************************************************************************/

void ofApp::setup()
{
    this->paused = true;
    this->advanceStep = false;
    
	ofSetLogLevel(OF_LOG_VERBOSE);
    ofSetVerticalSync(true);
    
    // This uses depth information for occlusion
    // rather than always drawing things on top of each other
    ofEnableDepthTest();
    
    // This sets the camera's distance from the object
    this->camera.setDistance(25);
    
    // Initialize from GL world:
    this->openCL.setupFromOpenGL();
    
    // Set the bounds of the simulation:
    AABB bounds((float3)(-5.0f, -5.0f, -5.0f), (float3)(5.0f, 5.0f, 5.0f));

    // Instantiate the simulator:
    this->simulation =
        std::shared_ptr<Simulation>(new Simulation(this->openCL, bounds));
}

void ofApp::update()
{
    if (this->isPaused()) {
        if (this->advanceStep) {
            this->simulation->step();
            this->advanceStep = false;
        }
    } else {
        this->simulation->step();
    }
}

void ofApp::draw()
{
    ofBackground(0);
    
    this->camera.begin();
 
        this->simulation->draw();

    this->camera.end();
    
    // Show the current frame rate and frame count
    ofSetColor(255);
    ofFill();
    ofDrawBitmapString(ofToString(ofGetFrameRate()) + " fps", 10, 15);
    ofDrawBitmapString("Frame: " + ofToString(this->simulation->getFrameNumber()), 10, 30);
    if (this->isPaused()) {
        ofDrawBitmapString("Paused", 10, 45);
    }
}

bool ofApp::isPaused() const
{
    return this->paused;
}

void ofApp::togglePaused()
{
    this->paused = !this->paused;
}

void ofApp::keyPressed(int key)
{
    switch (key) {
        // Pause
        case 'p':
        case ' ':
        {
            this->togglePaused();
        }
        break;
        // Step
        case 's':
        {
            this->advanceStep = true;
        }
        break;
        // Reset
        case 'r':
            break;
    }
}

void ofApp::keyReleased(int key)
{

}

void ofApp::mouseMoved(int x, int y )
{

}

void ofApp::mouseDragged(int x, int y, int button)
{

}

void ofApp::mousePressed(int x, int y, int button)
{

}

void ofApp::mouseReleased(int x, int y, int button)
{

}

void ofApp::windowResized(int w, int h)
{

}

void ofApp::gotMessage(ofMessage msg)
{

}

void ofApp::dragEvent(ofDragInfo dragInfo)
{

}

/******************************************************************************/
