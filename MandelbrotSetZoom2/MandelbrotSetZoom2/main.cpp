// press arrow keys to move the viewport 
//f1,f2 to change gradient of zoom in and out and to move, printed(printf) as diff 
//press 'w' to zoom in and press 'e' to zoom out
//press 'c' to change color 
//To zoom out u should increase gradient(diff) by pressing f2

#include <stdio.h>
#include <GL/freeglut.h>
#include <math.h>

extern "C" int cuComputeMandelbrotSet (int *ptr, int width, int height, GLdouble Rmin, GLdouble Rmax, GLdouble Imin, GLdouble Imax, int nIterations); 
extern "C" bool resetCUDADevice();

int nx, ny,c=0;
GLdouble diff = 0.1;
bool isCPU = true;
//GLdouble realMax=-1.49436,realMin=-1.4944f,imagMax=0.347133696469,imagMin=0.347130303531,realInc,imagInc;
GLdouble realMax=0.75f,realMin=-2.25f,imagMax=1.25f,imagMin=-1.25f,realInc,imagInc;
//GLdouble realMax = -0.774f, realMin = -0.782f, imagMax = 0.139f, imagMin = 0.133f, realInc, imagInc;
GLdouble r, g, b,tr,tg,tb;

typedef struct {
    GLdouble x, y;
} complex;

struct cuComplex 
{
    GLdouble   r;
    GLdouble   i;

    cuComplex( GLdouble a, GLdouble b ) : r(a), i(b)  {}

    GLdouble magnitude2( void ) { return r * r + i * i; }

    cuComplex operator*(const cuComplex& a) 
	{
        return cuComplex(r*a.r - i*a.i, i*a.r + r*a.i);
    }

    cuComplex operator+(const cuComplex& a) 
	{
        return cuComplex(r+a.r, i+a.i);
    }
};

/*complex complexSquare(complex c) {
    complex cSq;

    cSq.x = c.x * c.x - c.y * c.y;
    cSq.y = 2 * c.x * c.y;
    return (cSq);
}

int iterate (complex zInit,int maxIter)
{
    complex z = zInit;
    int cnt=0;

    while((z.x*z.x+z.y*z.y<=4.0)&&(cnt<maxIter)){
        z= complexSquare(z);
        z.x += zInit.x;
        z.y += zInit.y;
        cnt++;
    }
    return (cnt);
}
 */


int MandelbrotColor( int col, int row, int width, int height, GLdouble Rmin, GLdouble Rmax, GLdouble Imin, GLdouble Imax, int nIterations)
{
	// Normalize (col, row) to {(R,I) | Rmin < R < Rmax, Imin < I < Imax }
	GLdouble R = ((Rmax - Rmin)/(float)width) * (float)col + Rmin;
	GLdouble I = ((Imax - Imin)/(float)height) * (float)row + Imin;
    //float I = Imax - ((Imax - Imin)/(float)height) * (float)y;		// Bottom up

    cuComplex c(R, I);
    cuComplex a(R, I) ;
    int i = 0;
    for (i = 0; i < nIterations; i++) 
    {
        a = a * a + c;
        if (a.magnitude2() > 4.0)
            return i;
    }
   
    return i;
}

void ComputeMandelbrotSet(int width, int height, GLdouble Rmin, GLdouble Rmax, GLdouble Imax, GLdouble Imin, int maxIter)
{
	int cnt = 0; 
	for (int row = 0; row < height; row++) 
	{
        for (int col = 0; col < width; col++) 
		{
              cnt= MandelbrotColor (col, row, width, height, Rmin, Rmax, Imin, Imax, maxIter);

			  
            if (cnt == maxIter)
                {}
            else {                 
                       if(cnt>=0&&cnt<=31)   {b=cnt*4; g=cnt*8; r=0;  }
                  else if(cnt>=32&&cnt<=63)  {b=200; g=500-cnt*8; r=0;  }                      
		  else if(cnt>=64&&cnt<=95)  {b=200; g=0; r=(cnt-64)*4;}
		  else if(cnt>=96&&cnt<=127) {r=200; g=0; b=1000-cnt*8;}
		  else if(cnt>=128&&cnt<=159){r=200; g=(cnt-128)*8; b=0;}
		  else if(cnt>=160&&cnt<=191){g=200; r=1500-cnt*8; b=0;}
		  else if(cnt>=192&&cnt<=223){g=200; r=0; b=(cnt-192)*8;}
          else if(cnt>=224&&cnt<=255){g=230; r=(cnt-224)*8; b=256;}


                 //to change color by prssing key 'c'
                   	
                                tr=r;tb=b;tg=g;
                                switch(c)
				{
				case 0: break;
				case 1: r=tb;b=tr;break;
				case 2: r=tg;g=tr;break;
				case 3: b=tg;g=tb;break;
				case 4: r=tg; g=tb; b=tr; break;
				case 5: r=tb; g=tr; b=tg; break;
                      
				}	
			
			   
                    
                   glColor3f(r/256,g/256,b/256);
                   
                glVertex3d(col - nx / 2, row - ny / 2, 0.0f);
            }


			  //return index;
        }
    }
}

// Called to draw scene

void RenderScene(void) {
printf("Rendreing Scene == cpu? %d == %d == %d \n",isCPU,nx,ny);

    int maxIter;
                 r=g=b=256.0;    
    // Clear the window with current clearing color
    glClear(GL_COLOR_BUFFER_BIT);
    maxIter = 500;

    // Save matrix state and do the rotation
    glPushMatrix();

 //   realInc = (realMax - realMin) / (GLdouble)nx;
 //   imagInc = (imagMax - imagMin) / (GLdouble)ny;
    // Call only once for all remaining points
    glBegin(GL_POINTS);
    
	int *h_src = NULL;
	
	if (h_src != NULL) free (h_src);
		h_src = (int *) calloc (nx * ny, sizeof(int));
		
	

	if (h_src == NULL)
		h_src = (int *) calloc (nx * ny, sizeof(int));

	
	if(isCPU){
	 ComputeMandelbrotSet(nx, ny, realMin, realMax, imagMin, imagMax, maxIter);
	}else{
		resetCUDADevice();
		cuComputeMandelbrotSet(h_src, nx, ny, realMin, realMax, imagMin, imagMax, maxIter);

	int cnt = 0; 
	for (int row = 0; row < ny; row++) 
	{
        for (int col = 0; col < nx; col++) 
		{

			int index = (col + (ny-row-1) * nx);
			cnt = h_src[index];
			if (cnt == maxIter)
                {}
            else {                 
                       if(cnt>=0&&cnt<=31)   {b=cnt*4; g=cnt*8; r=0;  }
                  else if(cnt>=32&&cnt<=63)  {b=200; g=500-cnt*8; r=0;  }                      
		  else if(cnt>=64&&cnt<=95)  {b=200; g=0; r=(cnt-64)*4;}
		  else if(cnt>=96&&cnt<=127) {r=200; g=0; b=1000-cnt*8;}
		  else if(cnt>=128&&cnt<=159){r=200; g=(cnt-128)*8; b=0;}
		  else if(cnt>=160&&cnt<=191){g=200; r=1500-cnt*8; b=0;}
		  else if(cnt>=192&&cnt<=223){g=200; r=0; b=(cnt-192)*8;}
          else if(cnt>=224&&cnt<=255){g=230; r=(cnt-224)*8; b=256;}


                 //to change color by prssing key 'c'
                   	
                                tr=r;tb=b;tg=g;
                                switch(c)
				{
				case 0: break;
				case 1: r=tb;b=tr;break;
				case 2: r=tg;g=tr;break;
				case 3: b=tg;g=tb;break;
				case 4: r=tg; g=tb; b=tr; break;
				case 5: r=tb; g=tr; b=tg; break;
                      
				}	
			
			   
                    
                   glColor3f(r/256,g/256,b/256);
                   
                glVertex3d(col - nx / 2, row - ny / 2, 0.0f);
			}
		}
		}
	 free (h_src);
	}

	 /*

    for (x = 0; x < nx; x++) {
        for (y = 0; y < ny; y++) {


            //cnt=iterate(z,maxIter);
/*
         
            t2=z;
            cnt = 0;

            while ((t2.x * t2.x + t2.y * t2.y <=4) && (cnt < maxIter)) {
                t.x=t2.x*t2.x-t2.y*t2.y;
                t.y=2*t2.x*t2.y;
                t2 = t;
                t2.x +=z.x;
                t2.y +=z.y;
                cnt++;
            }	

	
			  cnt = MandelbrotColor (x, y, nx, ny, realMin, realMax, imagMin, imagMax, maxIter); 

			  /*
			
			// Normalize (col, row) to {(R,I) | Rmin < R < Rmax, Imin < I < Imax }
	float R = ((realMax - realMin)/(float)nx) * (float)x + realMin;
	float I = ((imagMax - imagMin)/(float)ny) * (float)y + imagMin;
    //float I = Imax - ((Imax - Imin)/(float)height) * (float)y;		// Bottom up

    cuComplex ci(R, I);
    cuComplex a(R, I) ;
    int i = 0;
    for (i = 0; i < maxIter; i++) 
    {
        a = a * a + ci;
        if (a.magnitude2() > 4)
            break;
    }
   
    cnt = i;
		
             
            if (cnt == maxIter)
                {}
            else {                 
                       if(cnt>=0&&cnt<=31)   {b=cnt*4; g=cnt*8; r=0;  }
                  else if(cnt>=32&&cnt<=63)  {b=200; g=500-cnt*8; r=0;  }                      
		  else if(cnt>=64&&cnt<=95)  {b=200; g=0; r=(cnt-64)*4;}
		  else if(cnt>=96&&cnt<=127) {r=200; g=0; b=1000-cnt*8;}
		  else if(cnt>=128&&cnt<=159){r=200; g=(cnt-128)*8; b=0;}
		  else if(cnt>=160&&cnt<=191){g=200; r=1500-cnt*8; b=0;}
		  else if(cnt>=192&&cnt<=223){g=200; r=0; b=(cnt-192)*8;}
          else if(cnt>=224&&cnt<=255){g=230; r=(cnt-224)*8; b=256;}


                 //to change color by prssing key 'c'
                   	
                                tr=r;tb=b;tg=g;
                                switch(c)
				{
				case 0: break;
				case 1: r=tb;b=tr;break;
				case 2: r=tg;g=tr;break;
				case 3: b=tg;g=tb;break;
				case 4: r=tg; g=tb; b=tr; break;
				case 5: r=tb; g=tr; b=tg; break;
                      
				}	
			
			   
                    
                   glColor3f(r/256,g/256,b/256);
                   
                glVertex3d(x - nx / 2, y - ny / 2, 0.0f);
            }
           
        



        }
    }
	*/
    printf("realMax=%.12lf\nrealMin=%.12lf\nimagMax=%.12lf\nimagMin=%.12lf\nimaginc=%.12lf\nrealinc=%.12lf\ndiff=%.12lf\n", realMax, realMin, imagMax, imagMin, imagInc, realInc, diff);


    // Done drawing points
    glEnd();

    // Restore transformations
    glPopMatrix();

    // Flush drawing commands
    glutSwapBuffers();
}

// This function does any needed initialization on the rendering
// context.

void SetupRC() {
    // Black background
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear (GL_COLOR_BUFFER_BIT);

    // Set drawing color to green
    glColor3f(0.0f, 1.0f, 0.0f);
}

void SpecialKeys(int key, int x, int y) {
    if (key == GLUT_KEY_UP){
        imagMax += diff;
        imagMin += diff;
    }
    if (key == GLUT_KEY_DOWN){
        imagMax -= diff;
        imagMin -= diff;
    }
        

    if (key == GLUT_KEY_RIGHT){
        realMax += diff;
        realMin += diff;
    }
    if (key == GLUT_KEY_LEFT){
        realMax -= diff;
        realMin -= diff;
    }
    if (key == GLUT_KEY_F1)
        diff /= 10.0f;

    if (key == GLUT_KEY_F2)
        diff *= 10.0f;


    // Refresh the Window
    glutPostRedisplay();
}

void proNormalKeys(unsigned char key, int x, int y) {
	
	GLdouble mul=0.833333333333333333333333;//this multiplier helps to zoom in x and y  properly 

	printf("at pronormal key");
	


	 if (key == 'g'||key == 'G'){
		 printf(" key pressed g");
		// getchar();
	
		 if (isCPU)
			 {
				 isCPU = false;
			 }
			
	}

	 if (key == 'c'||key == 'C'){
		 if (!isCPU)
			 {
				 isCPU = true;
			 }
			
	}
       
    if (key == 'w'||key == 'W'){
        if(realMax-realMin>3*diff){
            realMax -= diff;
            realMin += diff;
            imagMax -= diff*mul;
            imagMin += diff*mul;
        }
        else
            diff=diff/10;
    }


    if (key == 'e' || key == 'E'){
        realMax += diff;
        realMin -= diff;
        imagMax += diff*mul;
        imagMin -= diff*mul;
    }

    if (key == 'r' || key == 'R'){
	c++;
	if(c==6)c=0;
		

	 
    }
    // Refresh the Window
    glutPostRedisplay();
}
/*void ChangeSize(int w, int h)
	{
	// Calculate new clipping volume
	GLdouble windowWidth;
	GLdouble windowHeight;
          nx=300;
          ny=300;

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if(h == 0)
		h = 1;

	// Keep the square square
	if (w <= h) 
		{
		windowHeight = 200.0f*(GLdouble)h/(GLdouble)w;
		windowWidth = 200.0f;
		}
    else 
		{
		windowWidth = 200.0f*(GLdouble)w/(GLdouble)h;
		windowHeight = 200.0f;
		}

        // Set the viewport to be the entire window
        glViewport(0, 0, w, h);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

	// Set the clipping volume
	glOrtho(-200.0f, windowWidth, -200.0f, windowHeight, -200.0f, 200.0f);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

	}*/

void ChangeSize(int w, int h) {
    nx = w;
    ny = h;
    GLdouble nRange;

    // Prevent a divide by zero
    if (h == 0)
        h = 1;

    // Set Viewport to window dimensions
    glViewport(0, 0, w, h);

    // Reset projection matrix stack
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // Establish clipping volume (left, right, bottom, top, near, far)
    if (w <= h){
        nRange=nx/2;
        glOrtho(-nRange, nRange, -nRange * h / w, nRange * h / w, -nRange, nRange);
    }
    else{
        nRange=ny/2;
        glOrtho(-nRange * w / h, nRange * w / h, -nRange, nRange, -nRange, nRange);
    }
        
    // Reset Model view matrix stack
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
//	glutPostRedisplay();

}

int main(int argc, char* argv[]) {

	resetCUDADevice();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500, 500);
    glutCreateWindow("Mandelbrot Zoom");
    glutReshapeFunc(ChangeSize);
    glutSpecialFunc(SpecialKeys);
    glutKeyboardFunc(proNormalKeys);
    glutDisplayFunc(RenderScene);
    SetupRC();
    glutMainLoop();
	resetCUDADevice();

    return 0;
}
