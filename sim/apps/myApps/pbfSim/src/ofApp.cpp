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
    AABB bounds(float3(-5.0f, -10.0f, -5.0f), float3(5.0f, 10.0f, 5.0f));

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

/**
 * Draws a "heads up display" that shows the status of the simulation, as well
 * as some other pieces of pertinent information
 */
void ofApp::drawHeadsUpDisplay()
{
    // Show the current frame rate and frame count
    ofSetColor(255);
    ofFill();

    int textYOffset = 15;
    
    // FPS
    ofDrawBitmapString(ofToString(ofGetFrameRate()) + " fps", 10, textYOffset += 15);

    // Current frame
    ofDrawBitmapString("Frame: " + ofToString(this->simulation->getFrameNumber()), 10, textYOffset += 15);

    
    // Hotkeys:
    ofDrawBitmapString("Hotkeys: 's' = step, 'p' or space = toggle pause, 'r' = reset, 'g' = toggle grid, 'd' = toggle visual debugging"
                       , 10, textYOffset += 15);
    
    // Paused flag
    if (this->isPaused()) {
        ofDrawBitmapString("Paused", 10, textYOffset += 15);
    }
}

void ofApp::draw()
{
    ofBackground(0);
    
    this->camera.begin();
    this->simulation->draw();
    this->camera.end();
    
    this->drawHeadsUpDisplay();
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
            {
                this->simulation->reset();
            }
            break;
        // Toggle visual debugging:
        case 'd':
            {
                this->simulation->toggleVisualDebugging();
            }
            break;
        // Toggle drawing the spatial grid:
        case 'g':
            {
                this->simulation->toggleDrawGrid();
            }
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
