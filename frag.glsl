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
#define _x_ 0 /*'___' : macro name with a double underscore is reserved - unintented behavior is possible*/
#define RED 1
#define BLU 2
#define GRN 3
#define YLW 4
#define WHI 5
#define ORG 6
#define SHA 7
const vec3[] COLORS = {
	vec3(.03),
	vec3(1, 0, 0),
	vec3(0, 0, 1),
	vec3(0, 1, 0),
	vec3(1, 1, 0),
	vec3(1),
	vec3(1., .3, .0),
	vec3(.9, .9, .8),
};
int i;
const float ROUNDING = 2.5;
const float SIDE = 12. - ROUNDING;
const float UNIT = SIDE * 2. + 2. * ROUNDING;
const float SPACING = UNIT;
const vec3 off = vec3(0, SPACING, -SPACING);
const vec4 rot = vec4(0, HALFPI, -HALFPI, PI);
vec3 gHitPosition = vec3(0);
const vec3[] gCubeOff = {
	off.zyz, off.xyz, off.yyz, off.zxz, off.xxz, off.yxz, off.zzz, off.xzz, off.yzz,
	off.zyx, off.xyx, off.yyx, off.zxx, off.yxx, off.zzx, off.xzx, off.yzx, off.zyy,
	off.xyy, off.yyy, off.zxy, off.xxy, off.yxy, off.zzy, off.xzy, off.yzy,
};
const vec3[] gCubeRot = {
	rot.yxx, rot.xxx, rot.xxx, rot.yxx, rot.xxx, rot.zxx, rot.wxx, rot.wxx, rot.zxx,
	rot.xxz, rot.xyx, rot.zxz, rot.yyx, rot.zyx, rot.yxz, rot.xzx, rot.wxz, rot.yyx,
	rot.xyx, rot.xyx, rot.yyx, rot.xwx, rot.zyx, rot.wyx, rot.xwx, rot.xwx,
};
int gCubeCol1[26], gCubeCol2[26], gCubeCol3[26];
bool gCubeHidden[26];
int gHitIndex;
#define F 0 // front
#define G 1 // front(reverse)
#define L 2 // left
#define M 3 // left(reverse)
#define R 4 // right
#define S 5 // right(reverse)
#define U 6 // up
#define V 7 // up(reverse)
#define D 8 // down
#define E 9 // down(reverse)
#define B 10 // back
#define C 11 // back(reverse)
const int gNumMovements = 2;
const int[] gMovements = {
	F, F, V, D, L, M, F, F, M, S,
};
int gCurrentMovement;
float gCurrentMovementProgress;
const float MOVEMENT_TIME_SECONDS = 1.;

mat2 rot2(float a){float s=sin(a),c=cos(a);return mat2(c,s,-s,c);}

vec2 oneSidedCube(vec3 p, int cubeIndex)
{
	vec2 mc = vec2(length(max(abs(p) - SIDE, 0.)) - ROUNDING, 0);
	float tc = length(max(abs(p + vec3(0., 0., ROUNDING + .02)) - SIDE, 0.));

	return mc.x < tc ? mc : vec2(tc, float(gCubeCol1[cubeIndex]));
}

vec2 twoSidedCube(vec3 p, int cubeIndex)
{
	vec2 mc = oneSidedCube(p, cubeIndex);
	float fc = length(max(abs(p + vec3(0, -ROUNDING - .02, 0.)) - SIDE, 0.));

	return mc.x < fc ? mc : vec2(fc, float(gCubeCol2[cubeIndex]));
}

// also includes the shaft
vec2 centerCube(vec3 p, int cubeIndex, vec3 pos, vec3 rot)
{
	p -= pos;
	p.xy *= rot2(rot.x);
	p.yz *= rot2(rot.y);
	p.xz *= rot2(rot.z);
	vec2	mc = oneSidedCube(p, cubeIndex),
		sh = vec2(max(length(p.xy)-UNIT/4., length(max(abs(p-vec3(0.,0.,UNIT*.375)) - UNIT*.75, 0.))), 7);
	
	mc.x = max(mc.x, -(length(p+vec3(0.,0.,-UNIT))-UNIT*1.2));
	return mc.x < sh.x ? mc : sh;
}

vec2 middleCube(vec3 p, int cubeIndex, vec3 pos, vec3 rot)
{
	p -= pos;
	p.xy *= rot2(rot.x);
	p.yz *= rot2(rot.y);
	p.xz *= rot2(rot.z);
	vec2 res = twoSidedCube(p, cubeIndex);
	res.x = max(res.x, -(length(p+vec3(0., 62., -62.))-80.));
	res.x = min(res.x, length(max(abs(p - vec3(0., -7., 7.)) - vec3(7., 9., 9.), 0.)));
	return res;
}

vec2 cornerCube(vec3 p, int cubeIndex, vec3 pos, vec3 rot)
{
	p -= pos;
	p.xy *= rot2(rot.x);
	p.yz *= rot2(rot.y);
	p.xz *= rot2(rot.z);
	vec2 mc = twoSidedCube(p, cubeIndex);
	float sc = length(max(abs(p + vec3(-ROUNDING - .02, 0., 0.)) - SIDE, 0.));

	if (sc < mc.x) mc = vec2(sc, float(gCubeCol3[cubeIndex]));
	vec3 cubepos = p + vec3(9., 9., -9.);
	float bit = max(length(max(abs(cubepos) - vec3(7.), 0.)), length(p+vec3(UNIT,UNIT,-UNIT))-UNIT*1.18);
	mc.x = min(mc.x, bit);
	return mc;
}

vec2 map(vec3 p)
{
	vec2 res = vec2(9e9, 0);
	for (i = 0; i < 26; i++) {
		if (gCubeHidden[i]) {
			continue;
		}
		vec3 offset = gCubeOff[i];
		vec2 cub;
		vec3 pa = p;
		switch (gCurrentMovement) {
		case F:
			if (offset.y == off.y) {
				pa.xz *= rot2(HALFPI * gCurrentMovementProgress);
			}
			break;
		case V:
			break;
		}
		/*
		if (offset.x == off.z) {
			pa.yz *= rot2(iTime);
		}
		if (offset.x == off.y) {
			pa.yz *= rot2(-iTime);
		}
		if (offset.x == off.x) {
			pa.yz *= rot2(-iTime * 0.8);
		}
		if (offset.y == off.z) {
			//pa.xz *= rot2(iTime);
		}
		if (
			!(offset.x == 0. && offset.y == 0.) &&
			!(offset.x == 0. && offset.z == 0.) &&
			!(offset.y == 0. && offset.z == 0.)
		) {
			offset *= sin(iTime) + 2.;
		}
		*/
		vec3 rot = gCubeRot[i];
		if (gCubeCol2[i] == _x_) {
			cub = centerCube(pa, i, offset, rot);
		} else if (gCubeCol3[i] == _x_) {
			cub = middleCube(pa, i, offset, rot);
		} else {
			cub = cornerCube(pa, i, offset, rot);
		}
		if (cub.x < res.x) {
			res = cub;
			gHitIndex = i;
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
		gHitPosition = ro + rd * r.z;

		//p.y += 100.;
		//p.z -= 10.;
		//p.yz *= rot2(-.9);
		//p = mod(p, 30.) - 15.;
		//p.xy *= rot2(iTime/2.);
		//p.yz *= rot2(iTime/3.);
		//p.xz*=rot2(sin(p.z*0.2)*0.2+iTime);
		vec2 m = map(gHitPosition);
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

vec3 colorHit(vec4 result, vec3 rd)
{
	vec3 shade = vec3(0);
	int material = int(result.w);
	if (0 <= material && material <= 7) {
		shade = COLORS[material];
	}
	vec3 normal = norm(gHitPosition, result.y);
	// coloring magic from https://www.shadertoy.com/view/sdVczz
	float diffuse = max(0., dot(normal, -rd));
	float fresnel = pow(1. + dot(normal, rd), 4.);
	float specular = pow(max(dot(reflect(rd, normal), -rd), 0.), 30.);
	float ambientOcc = clamp(map(gHitPosition + normal * .05).x / .05, 0., 1.);
	float scat = smoothstep(0., 1., map(gHitPosition - rd * .4).x / .4); // "sub surface scattering"
	shade = mix(specular + shade * (ambientOcc + .2) * (diffuse + scat * .1), shade, fresnel);
	//shade = mix(background, shade, exp(-.002 * result.y * result.y * result.y));
	return shade;
}

#if shadertoy == 1
void mainImage(out vec4 c, in vec2 v)
#else
out vec4 c;
in vec2 v;
void main()
#endif
{
	gCubeCol1[ 0] = RED; gCubeCol2[ 0] = YLW; gCubeCol3[ 0] = BLU;
	gCubeCol1[ 0] = RED; gCubeCol2[ 0] = YLW; gCubeCol3[ 0] = BLU;
	gCubeCol1[ 1] = RED; gCubeCol2[ 1] = BLU; gCubeCol3[ 1] = _x_;
	gCubeCol1[ 2] = RED; gCubeCol2[ 2] = BLU; gCubeCol3[ 2] = WHI;
	gCubeCol1[ 3] = RED; gCubeCol2[ 3] = YLW; gCubeCol3[ 3] = _x_;
	gCubeCol1[ 4] = RED; gCubeCol2[ 4] = _x_; gCubeCol3[ 4] = _x_;
	gCubeCol1[ 5] = RED; gCubeCol2[ 5] = WHI; gCubeCol3[ 5] = _x_;
	gCubeCol1[ 6] = RED; gCubeCol2[ 6] = GRN; gCubeCol3[ 6] = YLW;
	gCubeCol1[ 7] = RED; gCubeCol2[ 7] = GRN; gCubeCol3[ 7] = _x_;
	gCubeCol1[ 8] = RED; gCubeCol2[ 8] = WHI; gCubeCol3[ 8] = GRN;
	gCubeCol1[ 9] = YLW; gCubeCol2[ 9] = BLU; gCubeCol3[ 9] = _x_;
	gCubeCol1[10] = BLU; gCubeCol2[10] = _x_; gCubeCol3[10] = _x_;
	gCubeCol1[11] = BLU; gCubeCol2[11] = WHI; gCubeCol3[11] = _x_;
	gCubeCol1[12] = YLW; gCubeCol2[12] = _x_; gCubeCol3[12] = _x_;
	gCubeCol1[13] = WHI; gCubeCol2[13] = _x_; gCubeCol3[13] = _x_;
	gCubeCol1[14] = GRN; gCubeCol2[14] = YLW; gCubeCol3[14] = _x_;
	gCubeCol1[15] = GRN; gCubeCol2[15] = _x_; gCubeCol3[15] = _x_;
	gCubeCol1[16] = WHI; gCubeCol2[16] = GRN; gCubeCol3[16] = _x_;
	gCubeCol1[17] = YLW; gCubeCol2[17] = ORG; gCubeCol3[17] = BLU;
	gCubeCol1[18] = BLU; gCubeCol2[18] = ORG; gCubeCol3[18] = _x_;
	gCubeCol1[19] = BLU; gCubeCol2[19] = ORG; gCubeCol3[19] = WHI;
	gCubeCol1[20] = YLW; gCubeCol2[20] = ORG; gCubeCol3[20] = _x_;
	gCubeCol1[21] = ORG; gCubeCol2[21] = _x_; gCubeCol3[21] = _x_;
	gCubeCol1[22] = WHI; gCubeCol2[22] = ORG; gCubeCol3[22] = _x_;
	gCubeCol1[23] = GRN; gCubeCol2[23] = ORG; gCubeCol3[23] = YLW;
	gCubeCol1[24] = ORG; gCubeCol2[24] = GRN; gCubeCol3[24] = _x_;
	gCubeCol1[25] = ORG; gCubeCol2[25] = GRN; gCubeCol3[25] = WHI;
	for (i = 0; i < 26; i++) {
		gCubeHidden[i] = false;
	}
	gCurrentMovement = -1;
	for (i = 0; i < gNumMovements; i++) {
		float until = float(i + 1) * MOVEMENT_TIME_SECONDS;
		if (float(until) < iTime) {
			int tmp;
			switch (gMovements[i]) {
			case F:
				tmp = gCubeCol1[2];
				gCubeCol1[2] = gCubeCol2[0];
				gCubeCol2[0] = gCubeCol2[17];
				gCubeCol2[17] = gCubeCol3[19];
				gCubeCol3[19] = tmp;
				tmp = gCubeCol2[11];
				gCubeCol2[11] = gCubeCol1[1];
				gCubeCol1[1] = gCubeCol1[9];
				gCubeCol1[9] = gCubeCol2[18];
				gCubeCol2[18] = tmp;
				tmp = gCubeCol2[19];
				gCubeCol2[19] = gCubeCol3[2];
				gCubeCol3[2] = gCubeCol1[0];
				gCubeCol1[0] = gCubeCol1[17];
				gCubeCol1[17] = tmp;
				break;
			case L:
				break;
			}
		} else {
			gCurrentMovement = gMovements[i];
			gCurrentMovementProgress = 1. - (float(until) - iTime) / MOVEMENT_TIME_SECONDS;
			break;
		}
	}

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
		vec3 shade = colorHit(result, rd);
		/*
		if (gHitIndex == 0) {
			gCubeHidden[gHitIndex] = true;
			result = march(gHitPosition, rd, 100); // TODO: how many steps?
			vec3 shade2 = result.x > 0. ? colorHit(result, rd) : col;
			shade = mix(shade, shade2, sin(iTime) * .5 + .5);
		}
		*/
		col = shade;
	}

	c = vec4(pow(col, vec3(.4545)), 1.0); // pow for gamma correction because all the cool kids do it
}
