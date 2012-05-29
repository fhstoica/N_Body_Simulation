#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <unistd.h>
#include <GL/glut.h>
#include <GL/glu.h>

const GLint dim = 3;
const GLint N   = 16;

GLfloat eye_x = 0;
GLfloat eye_y = 0;

typedef struct {
  GLfloat pos[N][dim];
  GLfloat vel[N][dim];
} ParticleCoords;

GLint globalFrameNo  = 0;
GLint globalFrameMax = 0;

void mouse_motion(GLint x, GLint y){
  eye_x = x;
  eye_y = y;
  glutPostRedisplay();
}

void Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ")
{
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while (std::string::npos != pos || std::string::npos != lastPos){
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void readParticleCoords(std::string f_in, std::vector<ParticleCoords*>& PosVec){
  
  std::ifstream fin;
  std::string c_line;
  std::vector<std::string> line_elements;
  GLint particle = 0;
  GLint col  = 0;
  GLint line = 0;
  fin.open(f_in.c_str());
  if (fin.is_open()){
    while( fin ){
      getline(fin, c_line);
      Tokenize(c_line, line_elements);      
      if( 2*dim*N < line_elements.size() ){
	ParticleCoords* const T = new ParticleCoords;
	col = 1;
	// We skip the first element. It is only the line number.
	std::vector<std::string>::iterator v_elem = line_elements.begin();
	line = atoi(v_elem->c_str());
	for(v_elem = line_elements.begin() + 1; v_elem != line_elements.end(); ++v_elem){	
	  if((col - 1)%(2*dim) == 0){
	    particle = static_cast<GLint>((col - 1)/(2*dim));
	  }
	  if((col - 1)%(2*dim) < dim){
	    T->pos[particle][(col-1)%(2*dim)] = atof(v_elem->c_str());
	  }
	  else{
	    T->vel[particle][(col-1)%(2*dim) - dim] = atof(v_elem->c_str());
	  }
	  ++col;
	}
	line_elements.clear();    
	PosVec.push_back(T);
      }
    }
  }
  fin.close();
}

void move(void){
    globalFrameNo += 1;
    globalFrameNo = globalFrameNo%globalFrameMax;
    glutPostRedisplay();
}

void create_particle_list(std::vector<ParticleCoords*>& PosVec, GLint FrameNo){
  
  glNewList(FrameNo, GL_COMPILE);
  glColor3f(1.0, 1.0, 0.0);
  glBegin(GL_POINTS);
      
  for(int partNo = 0; partNo < N; ++partNo){
    if(partNo%4 == 0){
      glColor3f(1.0, 1.0, 1.0);
    }
    else if(partNo%4 == 1){
      glColor3f(0.0, 1.0, 1.0);
    }
    else if(partNo%4 == 2){
      glColor3f(1.0, 0.0, 1.0);
    }
    else if(partNo%4 == 3){
      glColor3f(1.0, 1.0, 0.0);
    }    
    glVertex3f(PosVec[FrameNo]->pos[partNo][0], PosVec[FrameNo]->pos[partNo][1], PosVec[FrameNo]->pos[partNo][2]);      
  }
  glEnd();
  glColor3f(1.0, 1.0, 1.0);
  glEndList();
}

void display(void){
  glClear(GL_COLOR_BUFFER_BIT);/* Clears the window */
  glLoadIdentity();
  gluLookAt(10, 10, 25.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
  glPushMatrix();
  glRotatef(eye_x, 0, 1, 0);
  glRotatef(-eye_y, 0, 0, 1);
  glCallList(globalFrameNo);
  usleep(2000); 
  glPopMatrix();
  glutSwapBuffers();
  glFlush();
}

void keyboard(unsigned char key, int x, int y){
  /* Called when a key is pressed */
  if(key == 27){
    exit(0); /* 27 is the excape key */
  }
  else{
    printf("You pressed %c\n", key);
  }
}

void reshape(int w, int h){
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION); /* Select the projection matrix */
  glLoadIdentity(); /* Initialize it */
  gluPerspective(30, (GLfloat)w/(GLfloat)h, 0.0, 10.0);
  glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv){ 

  std::vector<ParticleCoords*> ParticlePositions;
  readParticleCoords( argv[1], ParticlePositions );

  globalFrameMax = ParticlePositions.size();

  glutInit(&argc, argv); /* Initializes OpenGL */
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
  glutInitWindowSize(800, 600); /* Set the window size */
  glutInitWindowPosition(100, 100); /* Set the window position */
  glutCreateWindow("Show_Cluster_Motion"); /* Create a window */

  for(unsigned int j = 1; j < ParticlePositions.size(); ++j){
    create_particle_list(ParticlePositions, j);
  }    
   
  glutDisplayFunc(display); /* Register the "display" function */
  glutKeyboardFunc(keyboard); /* Register the "keyboard" function */
  glutReshapeFunc(reshape); /* Register the "reshape" function */
  glutMotionFunc(mouse_motion); /* Register the "mouse_motion" function */
  glutIdleFunc(move); /* Register the "idle" function */
  glutMainLoop(); /* Enter the OpenGL main loop */
  return(0);
}

