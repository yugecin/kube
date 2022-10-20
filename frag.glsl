// vi: syntax=c
#version 430
#define iTime fpar[0].x
#define debugmov 1 //noexport
//#define flopineshade 1 //noexport
#define PI 3.14159265359
#define HALFPI 1.5707963268
#define fight(a,b) if(b.x<a.x)a=b; // union of 2 distance functions where x is distance and y is material id
layout (location=0) uniform vec4 fpar[2];
layout (location=2) uniform vec4 debug[2]; //noexport
layout (location=4) uniform sampler2D tex;
out vec4 c;
in vec2 v;
float g;
int i;

float _atan2(float y, float x)
{
	if(x>0.)return atan(y/x);
	if(x==0.)if(y>0.)return HALFPI;else return -HALFPI;
	if(y<0.)return atan(y/x)-PI;return atan(y/x)+PI;
}
float _atan2(vec2 v){return _atan2(v.y,v.x);}

float rand(vec2 p){return fract(sin(dot(p.xy,vec2(12.9898,78.233)))*43758.5453);}

mat2 rot2(float a){float s=sin(a),c=cos(a);return mat2(c,s,-s,c);}

float su(float d1, float d2, float k)
{
	float h = clamp(.5+.5*(d2-d1)/k,0.,1.);
	return mix(d2,d1,h)-k*h*(1.-h);
}
vec3 hitPosition = vec3(0);

// also includes the shaft
vec2 centerCube(vec3 p, int materialTop, vec3 pos, vec2 rot)
{
	p -= pos;
	vec2	mc = vec2(length(max(abs(p) - vec3(10.), 0.)) - 2, 0),
		tc = vec2(length(max(abs(p + vec3(0., 0., 2.02)) - vec3(10), 0.)), materialTop);
	
	fight(mc, tc);
	return mc;
}

vec2 twoSidedCube(vec3 p, int materialTop, int materialFront)
{
	vec2	mc = vec2(length(max(abs(p) - vec3(10.), 0.)) - 2, 0),
		tc = vec2(length(max(abs(p + vec3(0., 0., 2.02)) - vec3(10), 0.)), materialTop),
		fc = vec2(length(max(abs(p + vec3(0, -2.02, 0.)) - vec3(10), 0.)), materialFront);
	
	fight(mc, tc);
	fight(mc, fc);
	return mc;
}

vec2 middleCube(vec3 p, int materialTop, int materialFront, vec3 pos, vec2 rot)
{
	p -= pos;
	//p.xy *= rot2(rot.x);
	//p.xy *= rot2(rot.x);
	vec2 res = twoSidedCube(p, materialTop, materialFront);
	res.x = max(res.x, -(length(p+vec3(0., 62., -62.))-80.));
	res.x = min(res.x, length(max(abs(p - vec3(0., -7., 7.)) - vec3(7., 9., 9.), 0.)));
	return res;
}

vec2 cornerCube(vec3 p, int materialTop, int materialFront, int materialSide, vec3 pos, vec2 rot)
{
	p -= pos;
	rot *= HALFPI;
	p.xy *= rot2(rot.x);
	p.xy *= rot2(rot.y);
	vec2	mc = twoSidedCube(p, materialTop, materialFront),
		sc = vec2(length(max(abs(p + vec3(-2.02, 0., 0.)) - vec3(10), 0.)), materialSide);

	fight(mc, sc);
	vec3 cubepos = p + vec3(9., 9., -9.);
	float	bitcube = length(max(abs(cubepos) - vec3(7.), 0.)),
		bitpla1 = dot(cubepos, normalize(vec3(1, 0, -1))) - 4,
		bitpla2 = dot(cubepos, normalize(vec3(0, 1, -1))) - 4,
		bitpla3 = dot(cubepos, normalize(vec3(1, 1, 0))) - 4,
		bit = max(bitpla3, max(bitpla2, max(bitpla1, bitcube)));
	mc.x = min(mc.x, bit);
	return mc;
}

vec2 rotpos(vec2 res, vec3 pos, float b)
{
	return res;
}

float unit = 24. * 2.;
vec3 offs = vec3(0., unit, -unit);

vec2 cube(vec3 p)
{
	vec2	aaa = cornerCube(p, 1, 2, 3, offs.zyx, vec2(0., 1.)),
		aab = middleCube(p, 1, 2, offs.xyx, vec2(0.)),
		aac = cornerCube(p, 1, 2, 3, offs.yyx, vec2(0.)),
		bab = centerCube(p, 1, offs.xxx, vec2(0.));

	fight(aaa, aab);
	fight(aaa, aac);
	fight(aaa, bab);
	return aaa;
}

vec2 map(vec3 p)
{
	//vec2 c1 = middleCube(p, 1, 2);
	//vec2 c1 = cornerCube(p, 1, 2, 3);
	return cube(p);
	//return c1;
	//return vec2(length(p) - 5., 1);
}

vec3 norm(vec3 p, float dist_to_p)
{
	vec2 e = vec2(.001, 0);
	return normalize(dist_to_p - vec3(map(p-e.xyy).x, map(p-e.yxy).x, map(p-e.yyx).x));
}

// x=hit y=dist_to_p z=dist_to_ro w=material(if hit)
vec4 march(vec3 ro, vec3 rd, int maxSteps)
{
	vec4 r = vec4(0);
	for (i = 0; i < maxSteps && r.z < 350; i++){
		hitPosition = ro + rd * r.z;

		//p.y += 100.;
		//p.z -= 10.;
		//p.yz *= rot2(-.9);
		//p = mod(p, 30.) - 15.;
		//p.xy *= rot2(iTime/2.);
		//p.yz *= rot2(iTime/3.);
		//p.xz*=rot2(sin(p.z*0.2)*0.2+iTime);
		vec2 m = map(hitPosition);
		float distance = m.x;
		if (distance < .008) {
			r.x = 1;
			r.y = distance;
			r.w = m.y;
			break;
		}
		r.z += distance;
	}
	return r;
}

void main()
{
	vec2 q,lo,uv=v;uv.y/=1.77;

	vec3 ro = vec3(10 * sin(iTime), -30 * cos(iTime), -20);
	vec3 at = vec3(0, 0, -25);
	//vec3 ro = vec3(20, -50, -20);
	//ro.x = fpar[0].y/10.;
	//vec2 m = iMouse.xy/iResolution.xy;
	//ro.yz *= rot2(-m.y*PI+1.);
	//ro.xy *= rot2(-m.x*TAU);
	float horzAngle, vertAngle;
#if debugmov //noexport
	ro = debug[0].xyz; //noexport
	vertAngle = debug[1].y/20.; //noexport
	horzAngle = debug[1].x/20.; //noexport
#endif //noexport
	if (abs(vertAngle) < .001) {
		vertAngle = .001;
	}
	float xylen = sin(vertAngle);
	vertAngle = cos(vertAngle);
	at.x = ro.x + cos(horzAngle) * xylen;
	at.y = ro.y + sin(horzAngle) * xylen;
	at.z = ro.z + vertAngle;


	vec3 xx,l,cf=normalize(at-ro),
		cl=normalize(cross(cf,vec3(0,0,-1))),
		rd=mat3(cl,normalize(cross(cl,cf)),cf)*normalize(vec3(uv,1)),
		b;
	vec4 result = march(ro, rd, 200);

	vec3 col = vec3(.1) - length(uv) * .05;
	vec3 lightDirection = normalize(vec3(-.1, .4, -.3));
	lightDirection = -rd;

	if (result.x > 0) {
		// hit
		vec3 shade = vec3(0);
		float material = result.w;
		if (material == 0) {
			shade = vec3(.03);
		} else if (material == 1) {
			shade = vec3(1, 0, 0);
		} else if (material == 2) {
			shade = vec3(0, 0, 1);
		} else if (material == 3) {
			shade = vec3(0, 1, 0);
		}
		vec3 normal = norm(hitPosition, result.y);
		// coloring magic from https://www.shadertoy.com/view/sdVczz
		float diffuse = max(0., dot(normal, lightDirection));
		float fresnel = pow(1. + dot(normal, rd), 4.);
		float specular = pow(max(dot(reflect(-lightDirection, normal), -rd), 0.), 30.);
		float ambientOcc = clamp(map(hitPosition + normal * .05).x / .05, 0., 1.);
		float scat = smoothstep(0., 1., map(hitPosition + lightDirection * .4).x / .4);
		shade = mix(specular + shade * (ambientOcc + .2) * (diffuse + scat * .1), shade, fresnel);
		col = shade;
		//col = mix(col, shade, exp(-.002 * result.y * result.y * result.y));
	}

	c = vec4(pow(col, vec3(.4545)), 1.0); // pow for gamma correction because all the cool kids do it
}
