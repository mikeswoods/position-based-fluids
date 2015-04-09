/*******************************************************************************
 * ofApp
 * - The OpenFrameworks application entry point
 *
 * CIS563: Physcially Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#include "ofApp.h"
#include "Definitions.h"

/******************************************************************************/

using namespace std;

/******************************************************************************/

void ofApp::setup()
{
    ofSetVerticalSync(true);
    
    // This uses depth information for occlusion
    // rather than always drawing things on top of each other
    ofEnableDepthTest();
    
    // This sets the camera's distance from the object
    this->camera.setDistance(100);
    
    // Set the bounds of the simulation:
    AABB bounds(EigenVector3(-2.0f, 0.0f, -2.0f), EigenVector3(2.0f, 2.0f, 2.0f));
    
    // Instantiate the simulator:
    //this->simulation = std::make_shared<Simulation>(bounds, 100, 0.1f);
}

void ofApp::update()
{

}

void ofApp::draw()
{
    this->camera.begin();
    ofRotateX(ofRadToDeg(.5));
    ofRotateY(ofRadToDeg(-.5));
    
    ofBackground(0);
    
    ofSetColor(255,0,0);
    ofFill();
    ofDrawBox(30);
    ofNoFill();
    ofSetColor(0);
    ofDrawBox(30);
    
    ofPushMatrix();
    ofTranslate(0,0,20);
    ofSetColor(0,0,255);
    ofFill();
    ofDrawBox(5);
    ofNoFill();
    ofSetColor(0);
    ofDrawBox(5);
    ofPopMatrix();
    this->camera.end();
}

void ofApp::keyPressed(int key)
{

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
