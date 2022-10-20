
//#define dbg
#define messages // so it doesn't freeze when clicking (unnecessary but ok)
//#define nopopup // then screenshot works :^) (when also using "registerclass" and "messages")
//#define registerclass
//#define fpslimit
//#define nofullscreen
#ifndef XRES
#define XRES 1920
#endif
#ifndef YRES
#define YRES 1080
#endif

#define WIN32_LEAN_AND_MEAN
#define WIN32_EXTRA_LEAN
#include "windows.h"
#include <GL/gl.h>
#include <GL/glext.h>
char *vsh=
"#version 430\n"
"layout (location=0) in vec2 i;"
"out vec2 p;"
"out gl_PerVertex"
"{"
"vec4 gl_Position;"
"};"
"void main() {"
"gl_Position=vec4(p=i,0.,1.);"
"}"
;
#include "frag.glsl.c"

PIXELFORMATDESCRIPTOR pfd={0,1,PFD_SUPPORT_OPENGL|PFD_DOUBLEBUFFER, 32, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0};

#ifdef registerclass
WNDCLASSEX wc = {0};
#endif

/*
DEVMODE dmScreenSettings={ {0},0,0,sizeof(DEVMODE),0,DM_PELSWIDTH|DM_PELSHEIGHT,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1024,768,0,0,0,0,0,0,0,0,0,0};

/*
static DEVMODE dmScreenSettings = { {0},
    #if _MSC_VER < 1400
    0,0,148,0,0x001c0000,{0},0,0,0,0,0,0,0,0,0,{0},0,32,XRES,YRES,0,0,      // Visual C++ 6.0
    #else
    0,0,156,0,0x001c0000,{0},0,0,0,0,0,{0},0,32,XRES,YRES,{0}, 0,           // Visuatl Studio 2005
    #endif
    #if(WINVER >= 0x0400)
    0,0,0,0,0,0,
    #if (WINVER >= 0x0500) || (_WIN32_WINNT >= 0x0400)
    0,0
    #endif
    #endif
    };
    */

//#pragma code_seg(".main")
#ifdef nofullscreen
LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch (msg) {
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	default: return DefWindowProc(hwnd, msg, wParam, lParam);
	}
}
#endif
//gcc+ld? int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nShowCmd)
//gcc+link?  int WINAPI _WinMainCRTStartup(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nShowCmd)
//gcc+crinkler void mainCRTStartup(void)
//gcc+crinkler subsystem:windows
void WinMainCRTStartup(void)
{
#ifdef dbg
	char log[1024];
	int logsize;
#endif
#ifdef messages
	MSG msg;
#endif
#ifdef nopopup
	RECT rect;
#endif
	DEVMODE dm = {0};
	dm.dmSize = sizeof(DEVMODE);
	dm.dmFields = DM_PELSHEIGHT | DM_PELSWIDTH;
	dm.dmPelsWidth = XRES;
	dm.dmPelsHeight = YRES;

	float fparams[4*2];
	int it,t,
#ifdef fpslimit
		t2=-1004,
#endif
		k,tex;
#ifndef nofullscreen
	//ChangeDisplaySettings(&dm,CDS_FULLSCREEN);
#endif

#ifdef registerclass
	wc.cbSize = sizeof(WNDCLASSEX);
	wc.style = 0;
	wc.lpfnWndProc = DefWindowProc;
#ifdef nofullscreen
	wc.lpfnWndProc = WndProc;
#endif
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hInstance = (HINSTANCE)0x400000;
	wc.hIcon = LoadIcon(NULL, IDI_APPLICATION); /*large icon (alt tab)*/
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground = (HBRUSH) COLOR_WINDOW;
	wc.lpszMenuName = NULL;
	wc.lpszClassName = "metroclass";
	wc.hIconSm = LoadIcon(NULL, IDI_APPLICATION); /*small icon (taskbar)*/

	if (!RegisterClassEx(&wc)) {
		ExitProcess(0);
	}

	HANDLE hWnd = CreateWindowEx(WS_EX_APPWINDOW,wc.lpszClassName,"title",
#ifdef nofullscreen
		WS_OVERLAPPEDWINDOW |
#endif
#ifndef nopopup
		WS_POPUP |
#endif
		WS_VISIBLE, 0, 0, XRES, YRES, 0, 0, wc.hInstance, 0);
#else
	HANDLE hWnd = CreateWindow("static",0,
#ifndef nopopup
		WS_POPUP |
#endif
		WS_VISIBLE | WS_MAXIMIZE, 0, 0, XRES, YRES, 0, 0, 0, 0);
#endif
	HDC hDC = GetDC(hWnd);
	SetPixelFormat(hDC, ChoosePixelFormat(hDC, &pfd) , &pfd);
	wglMakeCurrent(hDC, wglCreateContext(hDC));
	/*
	GLuint p = ((PFNGLCREATEPROGRAMPROC)wglGetProcAddress("glCreateProgram"))();
	GLuint s = k =((PFNGLCREATESHADERPROC)(wglGetProcAddress("glCreateShader")))(GL_VERTEX_SHADER);
	((PFNGLSHADERSOURCEPROC)wglGetProcAddress("glShaderSource"))(s,1, &vsh,0);
	((PFNGLCOMPILESHADERPROC)wglGetProcAddress("glCompileShader"))(s);
	((PFNGLATTACHSHADERPROC)wglGetProcAddress("glAttachShader"))(p,s);
	s = ((PFNGLCREATESHADERPROC)
	wglGetProcAddress("glCreateShader"))(GL_FRAGMENT_SHADER);
	((PFNGLSHADERSOURCEPROC)wglGetProcAddress("glShaderSource"))(s,1, &fsh,0);
	((PFNGLCOMPILESHADERPROC)wglGetProcAddress("glCompileShader"))(s);
	((PFNGLATTACHSHADERPROC)wglGetProcAddress("glAttachShader"))(p,s);
	((PFNGLLINKPROGRAMPROC)wglGetProcAddress("glLinkProgram"))(p);
	((PFNGLUSEPROGRAMPROC)wglGetProcAddress("glUseProgram"))(p);
	*/
	GLuint p = ((PFNGLCREATESHADERPROGRAMVPROC)wglGetProcAddress("glCreateShaderProgramv"))(GL_VERTEX_SHADER, 1, &vsh);
	GLuint s = ((PFNGLCREATESHADERPROGRAMVPROC)wglGetProcAddress("glCreateShaderProgramv"))(GL_FRAGMENT_SHADER, 1, &fragSource);
	((PFNGLGENPROGRAMPIPELINESPROC)wglGetProcAddress("glGenProgramPipelines"))(1, &k);
	((PFNGLBINDPROGRAMPIPELINEPROC)wglGetProcAddress("glBindProgramPipeline"))(k);
	((PFNGLUSEPROGRAMSTAGESPROC)wglGetProcAddress("glUseProgramStages"))(k, GL_VERTEX_SHADER_BIT, p);
	((PFNGLUSEPROGRAMSTAGESPROC)wglGetProcAddress("glUseProgramStages"))(k, GL_FRAGMENT_SHADER_BIT, s);
#ifdef dbg
	logsize = 0;
	((PFNGLGETPROGRAMINFOLOGPROC)wglGetProcAddress("glGetProgramInfoLog"))(s, sizeof(log), &logsize, log);
	if (log[0] && logsize) {
		MessageBoxA(NULL, log, "hi", MB_OK);
		ExitProcess(1);
	}
#endif
#ifndef nofullscreen
	ShowCursor(0);
#endif

	it=GetTickCount();
	do
	{
		t=GetTickCount()-it;

#ifdef messages
		while (PeekMessage(&msg, 0, 0, 0, PM_REMOVE)) {
			if (msg.message == WM_QUIT) {
				goto done;
			}
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
#endif

#ifdef fpslimit
		if (t - t2 > 50) {
#endif

			// should technically also bind the texture ... but it works without ..
			fparams[0] = t/1000.0f;
			((PFNGLPROGRAMUNIFORM4FVPROC)wglGetProcAddress("glProgramUniform4fv"))(s, 0, 2, fparams);
			glRecti(1,1,-1,-1);
			SwapBuffers(hDC);
#ifdef fpslimit
			t2 = t;
		}
#endif

		// do your intro mainloop here
		// RenderIntro(MMTime.u.sample);

	} while (t < 30000 && !GetAsyncKeyState(VK_ESCAPE));
done:
#ifndef nofullscreen
	//ChangeDisplaySettings(0,0);
	ShowCursor(1);
#endif
	ExitProcess(0);
}
