#include <stdio.h>
#include <stdlib.h>
#include <assert.h>   

#include "communication_utils.h"    
#include "cellular_grid.h"  
#include "settings.h"

#ifdef SVG
#include <time.h>

FILE* svg = NULL;

/***************************** SVG saving functions *****************************/

void create_render(char* path_to_svg_folder, int width, int height){
    // Creating file name in the format "automata_<width>x<height>_<time>"
    assert(svg==NULL);
    char date_str[20];
    struct tm *date_and_time;

    time_t now = time(NULL);
    date_and_time = gmtime(&now);

    strftime(date_str, sizeof(date_str), "%Y-%m-%d_%H:%M:%S", date_and_time);

    char full_path[255];
    sprintf(full_path,"%s/automata_%dx%d_%s.svg",path_to_svg_folder,width,height,date_str);

    // Opening and initializing SVG
    svg = fopen(full_path,"w");
    fprintf(svg,"<svg version='1.1' viewBox='0 0 %d %d' xmlns='http://www.w3.org/2000/svg'>\n",width,height);
    fprintf(svg,"<rect width='%d' height='%d' x='0' y='0' fill='white'/>\n",width,height);
}

void render_generation(cell_point* points, int nb_points){
    for(int i=0; i<nb_points; i++){
        if(points[i].gen == 0)
            fprintf(svg,"<rect width='0' height='1' x='%d' y='%d' fill='black'><animate id='gen%d' attributeName='width' values='1' begin='0s;gen%d.end' dur='%s'/></rect>\n",points[i].x,points[i].y,points[i].gen,ITERATIONS-1,SVG_GEN_DURATION);
        else
            fprintf(svg,"<rect width='0' height='1' x='%d' y='%d' fill='black'><animate id='gen%d' attributeName='width' values='1' begin='gen%d.end' dur='%s'/></rect>\n",points[i].x,points[i].y,points[i].gen,points[i].gen - 1,SVG_GEN_DURATION);
    }
}

void finish_render(){
    fprintf(svg,"</svg>");
    fclose(svg);
    svg = NULL;
}

#endif

/***************************** X11 Display functions *****************************/

#ifdef X11
#include <X11/Xlib.h>
#include <unistd.h>   

#define NIL (0) 

struct X11_elements{
    Display *dpy;
    Window w;
    GC gc;
};

struct X11_elements* main_screen = NULL;

void create_render(char* unused_string, int width, int height){
    assert(main_screen==NULL);
    main_screen = malloc(sizeof(struct X11_elements));

    main_screen -> dpy = XOpenDisplay(NIL);
    assert(main_screen -> dpy);

    int blackColor = BlackPixel(main_screen -> dpy, DefaultScreen(main_screen -> dpy));
    int whiteColor = WhitePixel(main_screen -> dpy, DefaultScreen(main_screen -> dpy));

	main_screen -> w = XCreateSimpleWindow(main_screen -> dpy, DefaultRootWindow(main_screen -> dpy), 0, 0, width, height, 0, blackColor, blackColor);

	XSelectInput(main_screen -> dpy, main_screen -> w, StructureNotifyMask);

	XMapWindow(main_screen -> dpy, main_screen -> w);

	main_screen -> gc = XCreateGC(main_screen -> dpy, main_screen -> w, 0, NULL);

	XSetForeground(main_screen -> dpy, main_screen -> gc, whiteColor);

	for(;;) {
	    XEvent e;
	    XNextEvent(main_screen -> dpy, &e);
	    if (e.type == MapNotify)
		  break;
    }
}

void render_generation(cell_point* points, int nb_points){
    XClearWindow(main_screen -> dpy, main_screen -> w);
    int i;
    for (i=0; i<nb_points; i++){
        XDrawPoint(main_screen -> dpy, main_screen -> w, main_screen -> gc, points[i].x, points[i].y);
    }
    XFlush(main_screen -> dpy);
}

void finish_render(){
    XDestroyWindow(main_screen -> dpy, main_screen -> w);
    free(main_screen);
    main_screen = NULL;
}


#endif

#ifdef NORENDER

void create_render(char* path_to_svg_folder, int width, int height){}

void render_generation(cell_point* points, int nb_points){}

void finish_render(){}

#endif