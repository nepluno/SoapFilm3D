#include <iostream>
#include <sstream>

#ifdef __APPLE__
#include <GL/glew.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#ifdef WIN32
#include <Windows.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#endif

#include "Sim.h"
#include "MeshIO.h"

Sim g_sim(false);

struct SimControl
{
    bool headless;
    
    int win_w;
    int win_h;
    
    bool step;
    bool run;
    bool autoload;
    
    double view_theta;
    double view_alpha;
    double view_dist;
    
    int mouse_x;
    int mouse_y;
    
    bool ldrag;
    int ldrag_start_x;
    int ldrag_start_y;
    bool rdrag;
    int rdrag_start_x;
    int rdrag_start_y;
    Sim::RenderMode render_mode;
    
    int selection_mode;
    
} g_sc;

void renderBitmapString(float x, float y, float z, std::string s)
{
    glColor3f(0, 0, 0);
	glRasterPos3f(x, y, z);
	for (size_t i = 0; i < s.size(); i++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, s[i]);
}

void display()
{
    // object center and zoom
    Vec3d center(0, 0, 0);
    for (size_t i = 0; i < g_sim.vs()->mesh().nv(); i++)
        center += g_sim.vs()->pos(i);
    center /= g_sim.vs()->mesh().nv();
    Vec3d radius(0, 0, 0);
    for (size_t i = 0; i < g_sim.vs()->mesh().nv(); i++)
        for (size_t j = 0; j < 3; j++)
            radius[0] = std::max(radius[0], (g_sim.vs()->pos(i) - center)[0]);
    double min_d = std::max(std::max(radius[0], radius[1]), radius[2]) * 2.2;
    g_sc.view_dist = std::max(min_d, g_sc.view_dist);
    
    glClearColor(1, 1, 1, 1);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45, (double)g_sc.win_w / g_sc.win_h, 0.001, 50);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, -g_sc.view_dist, 0, 0, 0, 0, 0, 0, 1);
    glRotated(g_sc.view_alpha, 1, 0, 0);
    glRotated(g_sc.view_theta, 0, 0, 1);
    
    glBegin(GL_LINES);
    glColor3d(1, 0, 0);     glVertex3d(-2, 0, 0);    glVertex3d(2, 0, 0);
    glColor3d(0, 1, 0);     glVertex3d(0, -2, 0);    glVertex3d(0, 2, 0);
    glColor3d(0, 0, 1);     glVertex3d(0, 0, -2);    glVertex3d(0, 0, 2);
    glEnd();
    
    g_sim.render(g_sc.render_mode, Vec2d((double)g_sc.mouse_x / g_sc.win_w * 2 - 1, 1 - (double)g_sc.mouse_y / g_sc.win_h * 2), g_sc.selection_mode);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, g_sc.win_w, 0, g_sc.win_h, -1, 1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    std::stringstream ss;
    ss << "T = " << g_sim.time();
    std::string s = ss.str();
    renderBitmapString(20, g_sc.win_h - 20, 0, s);
    
    glutSwapBuffers();
}

void idle()
{
    if (g_sc.run || g_sc.step)
    {
        g_sc.step = false;
        if(g_sim.time() == 0.0)
            g_sim.stepOutput(g_sc.headless);
        g_sim.step();
        std::cout << "Finished step: T = " << g_sim.time() << std::endl;
        g_sim.stepOutput(g_sc.headless);
        if (g_sim.isFinished())
            exit(0);
        
        if (!g_sc.headless)
            glutPostRedisplay();
    }
    
    if (g_sc.autoload)
    {
        if (!g_sim.load(1))
            exit(0);
        g_sim.vs()->update_dbg_quantities();
        g_sim.stepOutput(g_sc.headless);
        
        if (!g_sc.headless)
            glutPostRedisplay();
    }
}

void keyboard(unsigned char k, int x, int y)
{
    if (k == 27 || k == 'q' || k == 'Q')
    {
        exit(0);
    } else if (k == ' ')
    {
        g_sc.run = !g_sc.run;
    } else if (k == 's' || k == 'S')
    {
        g_sc.step = true;
    } else if (k == 'm' || k == 'M')
    {
        g_sc.render_mode = (Sim::RenderMode)(((int)g_sc.render_mode + (k == 'm' ? 1 : -1)) % ((int)Sim::RM_COUNT));
        std::cout << "Render mode: " << (int)g_sc.render_mode << std::endl;
        glutPostRedisplay();
    } else if (k == 'v' || k == 'V')
    {
        g_sc.selection_mode = (k == 'v' ? (g_sc.selection_mode | Sim::SM_VERTEX) : (g_sc.selection_mode & ~Sim::SM_VERTEX));
        std::cout << "Mouse cursor selecting" << ((g_sc.selection_mode & Sim::SM_VERTEX) ? " vertices" : "") << ((g_sc.selection_mode & Sim::SM_EDGE) ? " edges" : "") << ((g_sc.selection_mode & Sim::SM_FACE) ? " faces" : "") << "." << std::endl;
        glutPostRedisplay();
    } else if (k == 'e' || k == 'E')
    {
        g_sc.selection_mode = (k == 'e' ? (g_sc.selection_mode | Sim::SM_EDGE) : (g_sc.selection_mode & ~Sim::SM_EDGE));
        std::cout << "Mouse cursor selecting" << ((g_sc.selection_mode & Sim::SM_VERTEX) ? " vertices" : "") << ((g_sc.selection_mode & Sim::SM_EDGE) ? " edges" : "") << ((g_sc.selection_mode & Sim::SM_FACE) ? " faces" : "") << "." << std::endl;
        glutPostRedisplay();
    } else if (k == 'f' || k == 'F')
    {
        g_sc.selection_mode = (k == 'f' ? (g_sc.selection_mode | Sim::SM_FACE) : (g_sc.selection_mode & ~Sim::SM_FACE));
        std::cout << "Mouse cursor selecting" << ((g_sc.selection_mode & Sim::SM_VERTEX) ? " vertices" : "") << ((g_sc.selection_mode & Sim::SM_EDGE) ? " edges" : "") << ((g_sc.selection_mode & Sim::SM_FACE) ? " faces" : "") << "." << std::endl;
        glutPostRedisplay();
    } else if (k == 'n' || k == 'N')
    {
        g_sim.showPrimitiveInfo();
    } else if (k == ']' || k == '}')
    {
        g_sim.load(k == ']' ? 1 : 10);
        g_sim.vs()->update_dbg_quantities();
        glutPostRedisplay();
    } else if (k == '[' || k == '{')
    {
        g_sim.load(k == '[' ? -1 : -10);
        g_sim.vs()->update_dbg_quantities();
        glutPostRedisplay();
    } else if (k == '.' || k == '>')
    {
        g_sim.load(k == '.' ? 100 : 1000);
        g_sim.vs()->update_dbg_quantities();
        glutPostRedisplay();
    } else if (k == ',' || k == '<')
    {
        g_sim.load(k == ',' ? -100 : -1000);
        g_sim.vs()->update_dbg_quantities();
        glutPostRedisplay();
    } else if (k == 'o' || k == 'O')
    {
        MeshIO::saveOBJ(*(g_sim.vs()), "mesh.obj");
    } else if (k == 'a' || k == 'A')
    {
        g_sc.autoload = true;
    }

}

void mouse(int b, int s, int x, int y)
{
    if (b == GLUT_LEFT_BUTTON && s == GLUT_DOWN)
    {
        g_sc.ldrag = true;
        g_sc.ldrag_start_x = x;
        g_sc.ldrag_start_y = y;
    } else if (b == GLUT_LEFT_BUTTON && s == GLUT_UP)
    {
        g_sc.ldrag = false;
    } else if (b == GLUT_RIGHT_BUTTON && s == GLUT_DOWN)
    {
        g_sc.rdrag = true;
        g_sc.rdrag_start_x = x;
        g_sc.rdrag_start_y = y;
    } else if (b == GLUT_RIGHT_BUTTON && s == GLUT_UP)
    {
        g_sc.rdrag = false;
    }
    
    glutPostRedisplay();
}

void motion(int x, int y)
{
    if (g_sc.ldrag)
    {
        g_sc.view_theta += (x - g_sc.ldrag_start_x) * 1.0;
        g_sc.view_alpha += (y - g_sc.ldrag_start_y) * 1.0;
        
        g_sc.ldrag_start_x = x;
        g_sc.ldrag_start_y = y;
    }
    if (g_sc.rdrag)
    {
        g_sc.view_dist *= pow(2.0, (y - g_sc.rdrag_start_y) * 0.01);
        
        g_sc.rdrag_start_x = x;
        g_sc.rdrag_start_y = y;
    }
    
    g_sc.mouse_x = x;
    g_sc.mouse_y = y;
    
    glutPostRedisplay();
}

void passiveMotion(int x, int y)
{
    g_sc.mouse_x = x;
    g_sc.mouse_y = y;
    
    glutPostRedisplay();
}

void reshape(int w, int h)
{
    g_sc.win_w = w;
    g_sc.win_h = h;
    
    glViewport(0, 0, w, h);
    
    glutPostRedisplay();
}

int main(int argc, char * argv[])
{
    if (argc != 2 && argc != 3)
    {
        std::cout << "Usage:\n\tSoapFilm3D option_file [save_outputs].\n" << std::endl;
        return 0;
    }
    
    // simulation setup
    g_sc.run = false;
    g_sc.step = false;
    g_sc.autoload = false;
    
    g_sc.win_w = 1600;
    g_sc.win_h = 1400;
    
    g_sc.view_theta = 0;
    g_sc.view_alpha = 0;
    g_sc.view_dist = 4;
    
    g_sc.mouse_x = 0;
    g_sc.mouse_y = 0;
    
    g_sc.ldrag = false;
    g_sc.ldrag_start_x = 0;
    g_sc.ldrag_start_y = 0;
    g_sc.rdrag = false;
    g_sc.rdrag_start_x = 0;
    g_sc.rdrag_start_y = 0;
    
    g_sc.render_mode = Sim::RM_TRANSPARENT;
    
    g_sc.selection_mode = 0;//Sim::SM_VERTEX | Sim::SM_EDGE | Sim::SM_FACE;
    
    g_sc.headless = (argc > 2 && strcmp(argv[2], "headless") == 0);

    if (!g_sc.headless)
    {
        // glut setup
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_MULTISAMPLE | GLUT_DEPTH);
        glutInitWindowPosition(0, 0);
        glutInitWindowSize(g_sc.win_w, g_sc.win_h);
        glutCreateWindow("Soap3D");

		if (GLEW_OK != glewInit()) {
			throw std::runtime_error("Failed to initialize GLEW\n");
		}
        
        glutKeyboardFunc(keyboard);
        glutMouseFunc(mouse);
        glutMotionFunc(motion);
        glutPassiveMotionFunc(passiveMotion);
        glutIdleFunc(idle);
        glutDisplayFunc(display);
        glutReshapeFunc(reshape);
        
        glEnable(GL_MULTISAMPLE_ARB);
        glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
        glEnable(GL_LINE_SMOOTH);
        glEnable(GL_POLYGON_SMOOTH);
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
        
    }
    
    bool success = g_sim.init(argv[1], argc > 2, g_sc.headless);
    if (!success)
        return 1;

    std::cout << "Initialization complete. Starting the simulation..." << std::endl;

    // main loop
    if (g_sc.headless)
    {
        g_sc.run = true;
        while (true)
            idle();
    } else
    {
        glutMainLoop();
    }
    
}

