// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

/**
 off-screen rendering on X11, using GL-X
 Using PBuffers, which are supported from GLX version 1.3 (OpenGL 1.2)
 
 https://en.wikipedia.org/wiki/Pixel_buffer
 */

#include <cstdio>
#include <X11/Xlib.h>
#include <GL/glx.h>

Display * dpy = 0;
GLXPbuffer pbuf = 0;
GLXContext glxContext;


unsigned OffScreen::openContext()
{
    dpy = XOpenDisplay(0);
    
    if ( dpy == 0 )
    {
        fprintf(stderr, "Unable to open a connection to the X server\n");
        return 1;
    }
    
    if ( ! glXQueryExtension(dpy, 0, 0) )
    {
        fprintf(stderr, "glX Extension not supported\n");
        return 2;
    }
    
#ifndef GLX_VERSION_1_3
    int minor, major;
    glXQueryVersion( dpy, &major, &minor );
    fprintf(stderr, "GLX version %i.%i detected\n", major, minor );
    fprintf(stderr, "GLX version 1.3 or better is required\n");
    return 3;
#endif
    return 0;
}


static GLXPbuffer makeBuffer(const int width, const int height, int)
{
    dpy = XOpenDisplay(0);
    
    if ( dpy == 0 )
    {
        fprintf(stderr, "Unable to open a connection to the X server\n");
        return 0;
    }
    
    if ( ! glXQueryExtension(dpy, 0, 0) )
    {
        fprintf(stderr, "glX Extension not supported\n");
        return 0;
    }
    
#ifndef GLX_VERSION_1_3
    int minor, major;
    glXQueryVersion( dpy, &major, &minor );
    fprintf(stderr, "GLX version %i.%i detected\n", major, minor );
    fprintf(stderr, "GLX version 1.3 or better is required\n");
    return 0;
#endif
    
    GLint screen = XDefaultScreen(dpy);
    
    GLint attribList[]= {
        GLX_RENDER_TYPE,   GLX_RGBA_BIT,
        GLX_DRAWABLE_TYPE, GLX_PBUFFER_BIT,
        GLX_RED_SIZE,      8,
        GLX_GREEN_SIZE,    8,
        GLX_BLUE_SIZE,     8,
        //GLX_DEPTH_SIZE,  16,
        None };
    
    /* config pointer for Frame Buffer */
    GLXFBConfig * FBConfig;
    /* number of FBConfigs returned */
    int FBConfig_Count;
    /* get frame buffer configuration */
    FBConfig = glXChooseFBConfig(dpy, screen, attribList, &FBConfig_Count);
    
    if ( FBConfig_Count == 0 || FBConfig == 0 )
    {
        fprintf(stderr, "glXChooseFBConfig returned NULL\n");
        return 0;
    }
    
    /* PBuffer Creation */
    int attribListP[]= {
        GLX_PRESERVED_CONTENTS, GL_TRUE,
        GLX_PBUFFER_WIDTH,  width,
        GLX_PBUFFER_HEIGHT, height,
        None };
    
    pbuf = glXCreatePbuffer(dpy, FBConfig[0], attribListP);
    
    glxContext = glXCreateNewContext(dpy, FBConfig[0], GLX_RGBA_TYPE, NULL, GL_TRUE);
    
    XFree(FBConfig);
    
    if ( glxContext == 0 )
    {
        fprintf(stderr, "glXCreateNewContext returned NULL\n");
        return 0;
    }
    
    if ( 0 == glXMakeCurrent(dpy, pbuf, glxContext) )
    {
        fprintf(stderr, "Cannot make the PBuffer current\n");
        return 0;
    }
    
    return pbuf;
}

unsigned OffScreen::openBuffer(const int width, const int height, int)
{
    GLXPbuffer buf = makeBuffer(width, height, 0);
    glViewport(0, 0, width, height);
    return buf;
}
    
void OffScreen::bindBuffer(unsigned id)
{
}

void OffScreen::releaseBuffer()
{
    glXDestroyPbuffer(dpy, pbuf);
}

void OffScreen::closeContext()
{
    releaseBuffer();
    glXDestroyContext(dpy, glxContext);
    XCloseDisplay(dpy);
}


