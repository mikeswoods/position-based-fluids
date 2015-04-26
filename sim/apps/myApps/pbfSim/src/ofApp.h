#pragma once

#include "ofMain.h"
#include "ofxGui.h"
#include "MSAOpenCL.h"
#include "Simulation.h"

/*******************************************************************************
 * OpenFrameworks base class
 ******************************************************************************/
 
class ofApp : public ofBaseApp
{
    private:
        void setSineAnim();
        void setRampAnim();
        void setCompressAnim();
        void doResetBounds();
        void setAnimPeriod(float& period);
        void setAnimAmp(float& amp);
    
    protected:
        // UI controls:
        ofxGuiGroup gui;
        ofxToggle toggleAnimateBounds;
        ofxToggle toggleAnimateBothSides;
        ofxButton selectSineAnim;
        ofxButton selectRampAnim;
        ofxButton selectCompressAnim;
        ofxButton resetBounds;
        ofxFloatSlider periodAnimSlider;
        ofxFloatSlider ampAnimSlider;

        // Flags
        bool paused;
        bool advanceStep;
    
        // Simulation, etc.
        ofEasyCam camera;
        msa::OpenCL openCL;
        std::shared_ptr<Simulation> simulation;
    
        void drawHeadsUpDisplay(ofEasyCam& camera);
    
	public:
		void setup();
		void update();
		void draw();
    
        bool isPaused() const;
        void togglePaused();

		void keyPressed(int key);
};
