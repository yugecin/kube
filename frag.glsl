// vi: syntax=c
#version 430
#define shadertoy 0
#if shadertoy == 0
#define iTime fpar[0].x
#define debugmov 1 //noexport
layout (location=0) uniform vec4 fpar[2];
layout (location=2) uniform vec4 debug[2]; //noexport
#endif
#define PI 3.14159265359
#define HALFPI 1.5707963268
#define fight(a,b) if(b.x<a.x)a=b; // union of 2 distance function resuls where x is distance and y is material id
#define _x_ 0. /*'___' : macro name with a double underscore is reserved - unintented behavior is possible*/
#define RED 1.
#define BLU 2.
#define GRN 3.
#define YLW 4.
#define WHI 5.
#define ORG 6.
#define SHA 7.
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

const float ROUNDING = 2.5;
const float SIDE = 12. - ROUNDING;
const float UNIT = SIDE * 2. + 2. * ROUNDING;
const float SPACING = UNIT;
const vec3 offs = vec3(0, SPACING, -SPACING);
const vec4 rot = vec4(0, HALFPI, -HALFPI, PI);

float su(float d1, float d2, float k)
{
	float h = clamp(.5+.5*(d2-d1)/k,0.,1.);
	return mix(d2,d1,h)-k*h*(1.-h);
}
vec3 hitPosition = vec3(0);

struct Cube {
	float col1, col2, col3;
	vec3 offset;
	vec3 rotation;
};
Cube cubes[26];

vec2 oneSidedCube(vec3 p, float materialTop)
{
	vec2	mc = vec2(length(max(abs(p) - SIDE, 0.)) - ROUNDING, 0),
		tc = vec2(length(max(abs(p + vec3(0., 0., ROUNDING + .02)) - SIDE, 0.)), materialTop);

	return mc.x < tc.x ? mc : tc;
}

vec2 twoSidedCube(vec3 p, float materialTop, float materialFront)
{
	vec2	mc = oneSidedCube(p, materialTop),
		fc = vec2(length(max(abs(p + vec3(0, -ROUNDING - .02, 0.)) - SIDE, 0.)), materialFront);

	return mc.x < fc.x ? mc : fc;
}

// also includes the shaft
vec2 centerCube(vec3 p, float materialTop, vec3 pos, vec3 rot)
{
	p -= pos;
	p.xy *= rot2(rot.x);
	p.yz *= rot2(rot.y);
	p.xz *= rot2(rot.z);
	vec2	mc = oneSidedCube(p, materialTop),
		sh = vec2(max(length(p.xy)-UNIT/4., length(max(abs(p-vec3(0.,0.,UNIT*.375)) - UNIT*.75, 0.))), 7);
	
	mc.x = max(mc.x, -(length(p+vec3(0.,0.,-UNIT))-UNIT*1.2));
	return mc.x < sh.x ? mc : sh;
}

vec2 middleCube(vec3 p, float materialTop, float materialFront, vec3 pos, vec3 rot)
{
	p -= pos;
	p.xy *= rot2(rot.x);
	p.yz *= rot2(rot.y);
	p.xz *= rot2(rot.z);
	vec2 res = twoSidedCube(p, materialTop, materialFront);
	res.x = max(res.x, -(length(p+vec3(0., 62., -62.))-80.));
	res.x = min(res.x, length(max(abs(p - vec3(0., -7., 7.)) - vec3(7., 9., 9.), 0.)));
	return res;
}

vec2 cornerCube(vec3 p, float materialTop, float materialFront, float materialSide, vec3 pos, vec3 rot)
{
	p -= pos;
	p.xy *= rot2(rot.x);
	p.yz *= rot2(rot.y);
	p.xz *= rot2(rot.z);
	vec2	mc = twoSidedCube(p, materialTop, materialFront),
		sc = vec2(length(max(abs(p + vec3(-ROUNDING - .02, 0., 0.)) - SIDE, 0.)), materialSide);

	fight(mc, sc);
	vec3 cubepos = p + vec3(9., 9., -9.);
	float bit = max(length(max(abs(cubepos) - vec3(7.), 0.)), length(p+vec3(UNIT,UNIT,-UNIT))-UNIT*1.18);
	mc.x = min(mc.x, bit);
	return mc;
}

vec2 map(vec3 p)
{
	vec2 res = vec2(9e9, 0);
	for (i = 0; i < 26; i++) {
		Cube c = cubes[i];
		vec2 cub;
		vec3 pa = p;
		vec3 of = c.offset;
		/*
		if (c.offset.x == offs.z) {
			pa.yz *= rot2(iTime);
		}
		if (c.offset.x == offs.y) {
			pa.yz *= rot2(-iTime);
		}
		if (c.offset.x == offs.x) {
			pa.yz *= rot2(-iTime * 0.8);
		}
		if (c.offset.y == offs.z) {
			//pa.xz *= rot2(iTime);
		}
		if (
			!(c.offset.x == 0 && c.offset.y == 0) &&
			!(c.offset.x == 0 && c.offset.z == 0) &&
			!(c.offset.y == 0 && c.offset.z == 0)
		) {
			of *= sin(iTime) + 2;
		}
		*/
		if (c.col2 == _x_) {
			cub = centerCube(pa, c.col1, of, c.rotation);
		} else if (c.col3 == _x_) {
			cub = middleCube(pa, c.col1, c.col2, of, c.rotation);
		} else {
			cub = cornerCube(pa, c.col1, c.col2, c.col3, of, c.rotation);
		}
		if (cub.x < res.x) {
			res = cub;
		}
	}
	return res;
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
	for (i = 0; i < maxSteps && r.z < 350.; i++){
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
		if (distance < .03) {
			r.x = 1.;
			r.y = distance;
			r.w = m.y;
			break;
		}
		r.z += distance;
	}
	return r;
}

#if shadertoy == 1
void mainImage(out vec4 c, in vec2 v)
#else
out vec4 c;
in vec2 v;
void main()
#endif
{
	cubes[ 0] = Cube(RED, YLW, BLU, offs.zyz, rot.yxx),
	cubes[ 1] = Cube(RED, BLU, _x_, offs.xyz, rot.xxx),
	cubes[ 2] = Cube(RED, BLU, WHI, offs.yyz, rot.xxx),
	cubes[ 3] = Cube(RED, YLW, _x_, offs.zxz, rot.yxx),
	cubes[ 4] = Cube(RED, _x_, _x_, offs.xxz, rot.xxx),
	cubes[ 5] = Cube(RED, WHI, _x_, offs.yxz, rot.zxx),
	cubes[ 6] = Cube(RED, GRN, YLW, offs.zzz, rot.wxx),
	cubes[ 7] = Cube(RED, GRN, _x_, offs.xzz, rot.wxx),
	cubes[ 8] = Cube(RED, WHI, GRN, offs.yzz, rot.zxx),
	cubes[ 9] = Cube(YLW, BLU, _x_, offs.zyx, rot.xxz),
	cubes[10] = Cube(BLU, _x_, _x_, offs.xyx, rot.xyx),
	cubes[11] = Cube(BLU, WHI, _x_, offs.yyx, rot.zxz),
	cubes[12] = Cube(YLW, _x_, _x_, offs.zxx, rot.yyx),
	cubes[13] = Cube(WHI, _x_, _x_, offs.yxx, rot.zyx),
	cubes[14] = Cube(GRN, YLW, _x_, offs.zzx, rot.yxz),
	cubes[15] = Cube(GRN, _x_, _x_, offs.xzx, rot.xzx),
	cubes[16] = Cube(WHI, GRN, _x_, offs.yzx, rot.wxz),
	cubes[17] = Cube(YLW, ORG, BLU, offs.zyy, rot.yyx),
	cubes[18] = Cube(BLU, ORG, _x_, offs.xyy, rot.xyx),
	cubes[19] = Cube(BLU, ORG, WHI, offs.yyy, rot.xyx),
	cubes[20] = Cube(YLW, ORG, _x_, offs.zxy, rot.yyx),
	cubes[21] = Cube(ORG, _x_, _x_, offs.xxy, rot.xwx),
	cubes[22] = Cube(WHI, ORG, _x_, offs.yxy, rot.zyx),
	cubes[23] = Cube(GRN, ORG, YLW, offs.zzy, rot.wyx),
	cubes[24] = Cube(ORG, GRN, _x_, offs.xzy, rot.xwx),
	cubes[25] = Cube(ORG, GRN, WHI, offs.yzy, rot.xwx);

	vec2 uv=v;uv.y/=1.77;

	vec3 ro;
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


	vec3	cf = normalize(at-ro),
		cl = normalize(cross(cf,vec3(0,0,-1))),
		rd = mat3(cl,normalize(cross(cl,cf)),cf)*normalize(vec3(uv,1)),
		col = vec3(.1) - length(uv) * .05;

	vec4 result = march(ro, rd, 300);

	if (result.x > 0.) { // hit
		vec3 shade = vec3(0);
		float material = result.w;
		if (material == _x_) {
			shade = vec3(.03);
		} else if (material == RED) {
			shade = vec3(1, 0, 0);
		} else if (material == BLU) {
			shade = vec3(0, 0, 1);
		} else if (material == GRN) {
			shade = vec3(0, 1, 0);
		} else if (material == YLW) {
			shade = vec3(1, 1, 0);
		} else if (material == WHI) {
			shade = vec3(1);
		} else if (material == ORG) {
			shade = vec3(1., .3, .0);
		} else if (material == SHA) {
			shade = vec3(.9, .9, .8);
		}
		vec3 normal = norm(hitPosition, result.y);
		// coloring magic from https://www.shadertoy.com/view/sdVczz
		float diffuse = max(0., dot(normal, -rd));
		float fresnel = pow(1. + dot(normal, rd), 4.);
		float specular = pow(max(dot(reflect(rd, normal), -rd), 0.), 30.);
		float ambientOcc = clamp(map(hitPosition + normal * .05).x / .05, 0., 1.);
		float scat = smoothstep(0., 1., map(hitPosition - rd * .4).x / .4); // "sub surface scattering"
		shade = mix(specular + shade * (ambientOcc + .2) * (diffuse + scat * .1), shade, fresnel);
		col = shade;
		//col = mix(col, shade, exp(-.002 * result.y * result.y * result.y));
	}

	c = vec4(pow(col, vec3(.4545)), 1.0); // pow for gamma correction because all the cool kids do it
}
