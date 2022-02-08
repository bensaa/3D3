//g++ main.cpp -o 3d3.bin `pkg-config gtkmm-3.0 --cflags --libs gio-unix-2.0` && ./3d3.bin
#include <cstdlib>
#include <iostream>
#include <gtkmm/drawingarea.h>
#include <gtkmm/application.h>
#include <gtkmm/window.h>
#include <gdkmm/pixbuf.h>
#include <glibmm/main.h>
#include <cairomm/context.h>
#include <giomm/resource.h>
#include <gdkmm/general.h> // set_source_pixbuf()
#include <glibmm/fileutils.h>
#include <math.h>
#include <vector>


Glib::RefPtr<Gtk::Application>* Happ = NULL;
float camerax = 1000,cameray = 100,cameraz = 1000;
float scrossvecx = 0,scrossvecy = 0,scrossvecz = -100;
float CRX = 0,CRY = 0;
float screenW,screenH;
float HscreenW,HscreenH;
float SRatio = 1.0;
float VAR = 1.4;
float VAL = -1.4;
float VAU = SRatio * 1.4;
float VAD = SRatio * -1.4;
float speedy = 0;

struct point{
    float sx = 0,sy = 0,sz = 0;
    float dx = 0,dy = 0, dz = 0;
    float fx = 0, fy = 0, fz = 0;
    float Rx = 0,Ry = 0;
    float TWAX = 0, TWAY = 0;
    float distance = 0;
};

struct kolmio{
    point* points[3];
    float r=255,g=0,b=0;
};

struct viisi{
    point* points[5];
    float r=255,g=0,b=0;
};

struct nelio{
    point* points[4];
    float r=255,g=0,b=0;
    float crossvecx = 0,crossvecy = 0,crossvecz = 0;
    float height = 0;
    float distsum = 0;
    int drawablesum = 0;
    bool draw = true;
};

struct esine{
    bool draw = true;
    std::vector<kolmio> kolmiot;
    std::vector<nelio> neliot;
    std::vector<point> points;
    float Rx = 0;
    float Ry = 0;
    float locx = 0,locy = 0,locz = 500;
};

std::vector<esine> esineet;
std::vector<nelio*> drawableNeliot;
std::vector<nelio*> goingtogetsplit;
std::vector<kolmio*> kolmiot;
std::vector<nelio*> neliot;
std::vector<viisi*> viitoset;
esine E;

nelio map[10000];
point mappoints[40000];

static void createmap(){
    for(int i = 0; i < 100;i++){
        for(int j = 0; j < 100;j++){
            map[i*100+j].g = 255;
            map[i*100+j].b = 0;
            map[i*100+j].r = 0;
            map[i*100+j].points[0] = &mappoints[i*400+j*4];
            map[i*100+j].points[0]->sx = j*100;
            map[i*100+j].points[0]->sz = (i+1)*100;
            map[i*100+j].points[1] = &mappoints[i*400+j*4 + 1];
            map[i*100+j].points[1]->sx = j*100;
            map[i*100+j].points[1]->sz = i*100;
            map[i*100+j].points[2] = &mappoints[i*400+j*4 + 2];
            map[i*100+j].points[2]->sx = (j+1)*100;
            map[i*100+j].points[2]->sz = i*100;
            map[i*100+j].points[3] = &mappoints[i*400+j*4 + 3];
            map[i*100+j].points[3]->sx = (j+1)*100;
            map[i*100+j].points[3]->sz = (i+1)*100;
        }
    }
	
    float heights[5]={200.0f,150.0f,100.0f,130.0f,180.0f};
    int coordx[5]= {25, 30, 75, 75, 50};
    int coordy[5]= {25, 30, 25, 75, 50};
    
    for(int i = 0;i < 5;i++){
        int steps = heights[i] /20;
        for(int i2 = 0; i2 < steps; i2++){
            for(int i3 = 0; i3 < steps; i3++){
                int index = (coordy[i]+i2)*100 + (coordx[i]+i3);
                if(i2 < i3){
                    map[index].height = (steps - i3)*20;   
                }
                else{
                    map[index].height = (steps - i2)*20;
                }
            }
        }
        for(int i2 = 0; i2 < steps; i2++){
            for(int i3 = 0; i3 < steps; i3++){
                int index = (coordy[i]-i2)*100 + (coordx[i]-i3);
                if(i2 < i3){
                    map[index].height = (steps - i3)*20;   
                }
                else{
                    map[index].height = (steps - i2)*20;
                }
            }
        }
    }
	
    for(int i=0;i<10000;i++){
        int v1 = map[i].height;
        int vn = map[i+100].height;
        int vne = map[i+101].height;
        int vnw = map[i+99].height;
        int vs = map[i-100].height;
        int vw = map[i-1].height;
        int ve = map[i+1].height;
        int vse = map[i-99].height;
        int vsw = map[i-101].height;
        map[i].points[0]->sy = (v1 + vn + vw + vnw)/4;
        map[i].points[1]->sy = (v1 + vw + vs + vsw)/4;
        map[i].points[2]->sy = (v1 + vs + ve + vse)/4;
        map[i].points[3]->sy = (v1 + vn + vne + ve)/4;
        if(map[i].height > 0){
            map[i].r = 255;
            map[i].g = 0;
        }
    }
    
}

static void dosplit(){
    for(int i = 0; i < goingtogetsplit.size();i++){
        nelio* n1 = goingtogetsplit[i];
        if(n1->drawablesum == 1){
            for(int i2 = 0; i2 < 4;i2++){
                if(n1->points[i2]->TWAX < VAR && n1->points[i2]->TWAX > VAL){
                    
                }
            }
        }
        else if(n1->drawablesum == 2){
        
        }
        else{
            
        }
    }
}

static void calcmap(){
    for(int j2 = 0;j2 < 10000;j2++){
        map[j2].distsum = 0;
        int visiblecorners = 0;
        for(int i = 0; i < 4;i++){
            point* P = map[j2].points[i];
            P->fz = P->sz - cameraz;
            float Tx = P->sx - camerax;
            float Ty = P->sy - cameray;
            float distance = sqrt(Tx * Tx + Ty * Ty + P->fz * P->fz );
			//Calculate Rotation Relative to camera position 
            P->distance = distance;
            if(distance < 2000){
                if(P->fz < 0 && Tx < 0){
                    P->TWAX = M_PI + asin((-1*Tx)/distance);
                    P->TWAX = P->TWAX - CRX;
                }
                else if(P->fz < 0){
                    P->TWAX = M_PI - asin(Tx/distance);
                    P->TWAX = P->TWAX - CRX;
                }
                else if(Tx < 0){
                    if(CRX < 0.5*M_PI){
                        P->TWAX = asin(Tx/distance);
                        P->TWAX = P->TWAX - CRX;
                    }
                    else{
                        P->TWAX = 2*M_PI + asin(Tx/distance);
                        P->TWAX = P->TWAX - CRX;
                    }
                }
                else if(CRX > 1.5 * M_PI){P->TWAX = 2*M_PI + asin(Tx/distance) - CRX;}

                else{P->TWAX = asin(Tx/distance) - CRX;}

                P->TWAY = asin(Ty/distance) - CRY;                        
				//Calculate Position on screen
                if(P->TWAX < VAR && P->TWAX > VAL){
                    visiblecorners++;
                    P->fx = HscreenW + (P->TWAX/VAR * HscreenW);
                    P->fy = HscreenH - (P->TWAY/VAR * HscreenW);
                    map[j2].distsum += distance;
                }
            }
        }
        if(visiblecorners == 0){
            map[j2].draw = false;
        }
        else if(visiblecorners == 4){
            map[j2].draw = true;
        }
        else{
            map[j2].draw = false;
            map[j2].drawablesum = visiblecorners;
            goingtogetsplit.emplace_back(&map[j2]);
        }
    }
    for(int i = 0; i < 10000;i++){
        nelio* KP = &(map[i]);
        float P1x = KP->points[1]->fx - KP->points[0]->fx; 
        float P1y = KP->points[1]->fy - KP->points[0]->fy; 
        float P1z = KP->points[1]->fz - KP->points[0]->fz; 
        float P2x = KP->points[3]->fx - KP->points[0]->fx; 
        float P2y = KP->points[3]->fy - KP->points[0]->fy; 
        float P2z = KP->points[3]->fz - KP->points[0]->fz; 

        KP->crossvecx = P1y*P2z - P1z*P2y;
        KP->crossvecy = P1z*P2x - P1x*P2z;
        KP->crossvecz = P1x*P2y - P1y*P2x;
        float value = scrossvecx * KP->crossvecx + scrossvecy * KP->crossvecy +scrossvecz * KP->crossvecz;

        if(value > 0 && KP->draw == true){
            drawableNeliot.emplace_back(KP);
        }
        else{
            KP->draw = false;
        }
    }
}
static void createCube(float x,float y,float z, float size){
    float hs = size/2;
    esineet.emplace_back();
    esineet[esineet.size()-1].locx = x;
    esineet[esineet.size()-1].locy = y;
    esineet[esineet.size()-1].locz = z;
//KULMAT
    esineet[esineet.size()-1].points.emplace_back((point){.sx = hs,.sy = hs,.sz = -hs});
    esineet[esineet.size()-1].points.emplace_back((point){.sx =-hs,.sy = hs,.sz = -hs});
    esineet[esineet.size()-1].points.emplace_back((point){.sx =-hs,.sy =-hs,.sz = -hs});
    esineet[esineet.size()-1].points.emplace_back((point){.sx = hs,.sy =-hs,.sz = -hs});
    esineet[esineet.size()-1].points.emplace_back((point){.sx = hs,.sy = hs,.sz =  hs});
    esineet[esineet.size()-1].points.emplace_back((point){.sx =-hs,.sy = hs,.sz =  hs});
    esineet[esineet.size()-1].points.emplace_back((point){.sx =-hs,.sy =-hs,.sz =  hs});
    esineet[esineet.size()-1].points.emplace_back((point){.sx = hs,.sy =-hs,.sz =  hs});
//SIVUT
    //Etu
    esineet[esineet.size()-1].neliot.emplace_back((nelio){&esineet[esineet.size()-1].points[0],&esineet[esineet.size()-1].points[1],
			&esineet[esineet.size()-1].points[2],&esineet[esineet.size()-1].points[3]});
    //Oikea
    esineet[esineet.size()-1].neliot.emplace_back((nelio){&esineet[esineet.size()-1].points[4],&esineet[esineet.size()-1].points[0],
			&esineet[esineet.size()-1].points[3],&esineet[esineet.size()-1].points[7]});
    //Lattia
    esineet[esineet.size()-1].neliot.emplace_back((nelio){&esineet[esineet.size()-1].points[7],&esineet[esineet.size()-1].points[3],
			&esineet[esineet.size()-1].points[2],&esineet[esineet.size()-1].points[6]});
    //Katto
    esineet[esineet.size()-1].neliot.emplace_back((nelio){&esineet[esineet.size()-1].points[0],&esineet[esineet.size()-1].points[4],
			&esineet[esineet.size()-1].points[5],&esineet[esineet.size()-1].points[1]});
    //Vasen
    esineet[esineet.size()-1].neliot.emplace_back((nelio){&esineet[esineet.size()-1].points[1],&esineet[esineet.size()-1].points[5],
			&esineet[esineet.size()-1].points[6],&esineet[esineet.size()-1].points[2]});
    //Taka
    esineet[esineet.size()-1].neliot.emplace_back((nelio){&esineet[esineet.size()-1].points[5],&esineet[esineet.size()-1].points[4],
			&esineet[esineet.size()-1].points[7],&esineet[esineet.size()-1].points[6]});
//Calculate X
    for(int j = 0; j < esineet[esineet.size()-1].points.size(); j++){
        point* P = &esineet[esineet.size()-1].points[j];
        float RXLength = sqrt(P->sx * P->sx + P->sz * P->sz);
        float Tsx = P->sx;
        float Tsy = P->sy;
        float Tsz = P->sz;
        if(P->sx < 0){Tsx *= -1;}
        if(P->sy < 0){Tsy *= -1;}
        if(P->sz < 0){Tsz *= -1;}
        if(P->sx < 0 && P->sz < 0){P->Rx = atan(Tsx/Tsz) + M_PI; }
        else if(P->sx < 0){P->Rx = 2*M_PI - atan(Tsx/Tsz);}
        else if(P->sz < 0){P->Rx = M_PI - atan(Tsx/Tsz);}
        else{P->Rx = atan(Tsx/Tsz);}
        float TRx = P->Rx;
        if(TRx > 2* M_PI){TRx = TRx - 2*M_PI;}
        P->dx = sin(TRx) * RXLength;
        P->dz = cos(TRx) * RXLength;
        if(P->dz < 0){P->dz *= -1;}
//Calculate Y
        Tsz = P->dz;
        if(P->dz < 0){Tsz *= -1;};
        float RYLength = sqrt(Tsz * Tsz + Tsy * Tsy);
        if(TRx > 0.5*M_PI && TRx < 1.5*M_PI){
            if(P->sy >= 0){P->Ry = M_PI - atan(Tsy/Tsz);}
            if(P->sy < 0){P->Ry = atan(Tsy/Tsz) + M_PI;}
        }
        else if(P->sy < 0){P->Ry = 2 * M_PI - atan(Tsy/Tsz);}
        else {P->Ry = atan(Tsy/Tsz);}

        float TRy = P->Ry;
        if(TRy > 2* M_PI){
            TRy = TRy - 2*M_PI;
        }
        P->dy = sin(TRy) * RYLength;
        P->dz = cos(TRy) * RYLength;
        P->fz = P->dz + esineet[j].locz - cameraz;
    }
}
static bool onKeyRelease(GdkEventKey* event){
    switch(event->hardware_keycode){
	case 38: {break;}
	case 25: {break;}
	case 40: {break;}
	case 39: {break;}
	case 65: {break;}
    }
    return true;
}
static void calcplace2(){
    drawableNeliot.clear();
    for(int j = 0;j < esineet.size();j++){
        for(int j2 = 0;j2 < esineet[j].neliot.size();j2++){
            esineet[j].neliot[j2].draw = true;
            esineet[j].neliot[j2].distsum = 0;
            for(int i = 0; i < 4;i++){
                point* P = esineet[j].neliot[j2].points[i];
                P->fz = P->dz + esineet[j].locz - cameraz;
                float Tx = (esineet[j].locx-camerax)+P->dx;
                float Ty = (esineet[j].locy-cameray)+P->dy;
                float distance = sqrt(Tx * Tx + Ty * Ty + P->fz * P->fz );
				//Calculate Rotation Relative to camera position 
                P->distance = distance;
                if(P->fz < 0 && Tx < 0){
                    P->TWAX = M_PI + asin((-1*Tx)/distance);
                    P->TWAX = P->TWAX - CRX;
                }
                else if(P->fz < 0){
                    P->TWAX = M_PI - asin(Tx/distance);
                    P->TWAX = P->TWAX - CRX;
                }
                else if(Tx < 0){
                    if(CRX < 0.5*M_PI){
                        P->TWAX = asin(Tx/distance);
                        P->TWAX = P->TWAX - CRX;
                    }
                    else{
                        P->TWAX = 2*M_PI + asin(Tx/distance);
                        P->TWAX = P->TWAX - CRX;
                    }
                }
                else if(CRX > 1.5 * M_PI){
                    P->TWAX = 2*M_PI + asin(Tx/distance) - CRX;
                }
                else{
                    P->TWAX = asin(Tx/distance) - CRX;
                }
                P->TWAY = asin(Ty/distance) - CRY;         

				//Calculate Position on screen
                if(P->TWAX > VAL && P->TWAX < VAR){
                    if(P->TWAY > VAD && P->TWAY < VAU){
                        esineet[j].draw = true;
                        P->fx = HscreenW + (P->TWAX/VAR * HscreenW);
                        P->fy = HscreenH - (P->TWAY/VAR * HscreenW);
                        esineet[j].neliot[j2].distsum += distance;
                    }
                    else{
                        esineet[j].draw = false;
                        break;
                    }
                }
                else{
                    esineet[j].draw = false;
                    break;
                }
            }
        }
    }
}

static void drawOrder(){
    for(int j = 0; j < drawableNeliot.size(); j++){
        for(int i = 0; i < drawableNeliot.size()-1; i++){
            if(drawableNeliot[i]->distsum < drawableNeliot[i+1]->distsum){
                nelio* temp = drawableNeliot[i];
                drawableNeliot[i] = drawableNeliot[i+1];
                drawableNeliot[i+1] = temp;
            }
        }
    }
}
static void checkVisibility(){
    for(int j = 0;j < esineet.size();j++){
        for(int i = 0; i < esineet[j].neliot.size();i++){
            nelio* KP = &(esineet[j].neliot[i]);
            float P1x = KP->points[1]->fx - KP->points[0]->fx; 
            float P1y = KP->points[1]->fy - KP->points[0]->fy; 
            float P1z = KP->points[1]->fz - KP->points[0]->fz; 
            float P2x = KP->points[3]->fx - KP->points[0]->fx; 
            float P2y = KP->points[3]->fy - KP->points[0]->fy; 
            float P2z = KP->points[3]->fz - KP->points[0]->fz; 
            
            KP->crossvecx = P1y*P2z - P1z*P2y;
            KP->crossvecy = P1z*P2x - P1x*P2z;
            KP->crossvecz = P1x*P2y - P1y*P2x;
            float value = scrossvecx * KP->crossvecx + scrossvecy * KP->crossvecy +scrossvecz * KP->crossvecz;
            
            if(value > 0 && KP->draw == true){
                drawableNeliot.emplace_back(&esineet[j].neliot[i]);
            }
            else{
                KP->draw = false;
            }
        }
    }
}

class Window3d : public Gtk::DrawingArea{
public:
    Window3d(){
        Glib::signal_timeout().connect( sigc::mem_fun(*this, &Window3d::on_timeout), 100 );
    };
    bool on_draw(const Cairo::RefPtr<Cairo::Context>& cr) override{
        calcplace2();
        drawOrder();
        calcmap();
        checkVisibility();
        cr->set_line_width(2.0);
        for(int i = 0; i < drawableNeliot.size();i++){
            nelio* KP = (drawableNeliot[i]);
            if(KP->draw){
                cr->set_source_rgba(KP->r, KP->g, KP->b, 1);
                cr->move_to(KP->points[0]->fx, KP->points[0]->fy);
                cr->line_to(KP->points[1]->fx, KP->points[1]->fy);
                cr->line_to(KP->points[2]->fx, KP->points[2]->fy);
                cr->line_to(KP->points[3]->fx, KP->points[3]->fy);
                cr->close_path();
                cr->fill_preserve();
                cr->set_source_rgb(0.0, 0.0, 0.0);
                cr->stroke();
            }
        }

        return true;
        
    };
    bool on_timeout(){
        cameray = cameray + speedy;
        if(cameray > 100.0f){
            speedy = speedy -5.0f;
        }
        else{
            cameray = 100.0f;
            speedy = 0.0f;
        }
        this->queue_draw();
        return true;
    };
};
Window3d* Hwin = NULL;

static bool onKeyPress(GdkEventKey* event){
    switch(event->hardware_keycode){
	case 38:{
		CRX = CRX - 0.05;
		if(CRX < 0){CRX = CRX + 2*M_PI;}
		break;
	}
	case 25:{
		float TCAX = sin(CRX);
		float TCAY = cos(CRX);
		camerax += TCAX * 30;
		cameraz += TCAY * 30;
		break;
	}
	case 39: {
		float TCAX = sin(CRX);
		float TCAY = cos(CRX);
		camerax -= TCAX * 30;
		cameraz -= TCAY * 30;
		break;
	}
	case 40:{
		CRX = CRX + 0.05;
		if(CRX >= 2*M_PI){CRX = CRX - 2*M_PI;}
		break;
	}
        
	case 31: {
		CRY = CRY + 0.05;
		if(CRY >= 0.5*M_PI){CRY = 0.5*M_PI;}
		break;
	}
	case 44:{
		E.Ry = E.Ry + 0.05;
		if(E.Ry > 2* M_PI){E.Ry = E.Ry - 2*M_PI;}
		break;
	}
	case 45:{
		CRY = CRY - 0.05;
		if(CRY <= -0.5*M_PI){CRY = -0.5*M_PI;}
		break;
	}
	case 46:{
		E.Ry = E.Ry - 0.05;
		if(E.Ry < 0){E.Ry = E.Ry + 2*M_PI;}
		break;
	}
	case 65: {
		speedy = speedy + 75;
		break;
	}
	case 56: {
	}
	case 9:{
		((Gtk::Application*)Happ)->quit();
	}
    }
    Hwin->queue_draw();
    return true;
}
static void onResize(){
    screenW = Hwin->get_width();
    screenH = Hwin->get_height();
    HscreenW = screenW/2;
    HscreenH = screenH/2;
    SRatio = screenH/screenW;
}

int main(int argc, char** argv){
    
    std::cout << asin(-100.0f/120.0f) << std::endl;
    createmap();
    createCube(1300,100,5000, 200);
    createCube(1000,100,5000, 200);
    createCube(700,100,5000, 200);

    auto app = Gtk::Application::create(argc, argv, "org.gtkmm.example");
    Happ = &app;
    Gtk::Window win;
    win.set_default_size(600,600);
    win.set_title("Window3d");
    Window3d drawarea_3D;
    Hwin = &drawarea_3D;
    win.fullscreen();
    win.signal_key_press_event().connect(sigc::ptr_fun(onKeyPress),false );
    win.signal_key_release_event().connect(sigc::ptr_fun(onKeyRelease),false );
    win.signal_check_resize().connect(sigc::ptr_fun(onResize),false );
    win.add(drawarea_3D);

    drawarea_3D.show();

    return app->run(win);
}

