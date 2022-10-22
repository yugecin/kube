// vi: syntax=c
#version 430
#define shadertoy 0
#if shadertoy == 0
#define iTime fpar[0].x
#define debugmov 0 //noexport
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
float gRounding, gSide, gUnit, gOffsetStuff;
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
int gCubeCol[3][26];
float gCubeOpacity[26];
bool gCubeHidden[26];
int gHitIndex = 0;
int gExclusiveCube = -1; // if set, only that cube will be considered.
                         //useful for consecutive map() calls for lighting after we already know what cube gets hit
#define F 0 // front
#define L 2 // left
#define R 4 // right
#define U 6 // up
#define D 8 // down
#define B 10 // back
const int gNumMovements = 12;
const int[] gMovements = {
	F, F, B, B, L, L, R, R, U, U, D, D
};
int gCurrentMovement;
float gCurrentMovementProgress;
const float MOVEMENT_TIME_SECONDS = .3;
const float HIDE_TIME_SECONDS = .2;
const float FADE_TIME_SECONDS = .4;
const int[] gCubeHiddenOrder = {
    9,
    5,
    25,
    24,
    17,
    11,
    3,
    6,
    8,
    2,
    14,
    19,
    1,
    20,
    22,
    16,
    23,
    7,
    18,
    //
    4,
    10,
    12,
    13,
    15,
    21,
    //0,
};
bool gHackFadeStuff;
bool gShaft;
float gTimeMod;

mat2 rot2(float a){float s=sin(a),c=cos(a);return mat2(c,s,-s,c);}

vec2 oneSidedCube(vec3 p, int cubeIndex)
{
	vec2 mc = vec2(length(max(abs(p) - gSide, 0.)) - gRounding, 0);
	float tc = length(max(abs(p + vec3(0., 0., gRounding + .02)) - gSide, 0.));

	return mc.x < tc ? mc : vec2(tc, float(gCubeCol[0][cubeIndex]));
}

vec2 twoSidedCube(vec3 p, int cubeIndex)
{
	vec2 mc = oneSidedCube(p, cubeIndex);
	float fc = length(max(abs(p + vec3(0, -gRounding - .02, 0.)) - gSide, 0.));

	return mc.x < fc ? mc : vec2(fc, float(gCubeCol[1][cubeIndex]));
}

// also includes the shaft
vec2 centerCube(vec3 p, int cubeIndex, vec3 pos, vec3 rot)
{
	p -= pos;
	p.xy *= rot2(rot.x);
	p.yz *= rot2(rot.y);
	p.xz *= rot2(rot.z);
	vec2	mc = oneSidedCube(p, cubeIndex),
		sh = vec2(max(length(p.xy)-gUnit/4., length(max(abs(p-vec3(0.,0.,gUnit*.375)) - gUnit*.75, 0.))), 7);
	
	mc.x = max(mc.x, -(length(p+vec3(0.,0.,-gUnit))-gUnit*1.2));
	return mc.x < sh.x || !gShaft ? mc : sh;
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
	float sc = length(max(abs(p + vec3(-gRounding - .02, 0., 0.)) - gSide, 0.));

	if (sc < mc.x) mc = vec2(sc, float(gCubeCol[2][cubeIndex]));
	vec3 cubepos = p + vec3(9., 9., -9.);
	float bit = max(length(max(abs(cubepos) - vec3(7.), 0.)), length(p+vec3(gUnit,gUnit,-gUnit))-gUnit*1.18);
	mc.x = min(mc.x, bit);
	return mc;
}

vec2 map(vec3 p)
{
	float tt = clamp(gTimeMod / 9., 0., 1.) * PI * 2;
	p.xy *= rot2(tt * .9 + PI * 2. * .1);
	if (gTimeMod >= 9.) {
		p.xy *= rot2(PI * 2. * .1 * clamp(gTimeMod - 9., 0., 1.));
	}
	//p.zy *= rot2(tt);
	float boundingbox = length(max(abs(p) - vec3(gUnit * 2.1), 0.));
	if (boundingbox > .2) {
		return vec2(boundingbox, 0.);
	}
	vec2 res = vec2(9e9, 0.);
	for (i = 0; i < 26; i++) {
		if (gCubeHidden[i] || (gExclusiveCube != -1 && gExclusiveCube != i)) {
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
		case L:
			if (offset.x == off.z) {
				pa.yz *= rot2(HALFPI * gCurrentMovementProgress);
			}
			break;
		case R:
			if (offset.x == off.y) {
				pa.zy *= rot2(HALFPI * gCurrentMovementProgress);
			}
			break;
		case B:
			if (offset.y == off.z) {
				pa.zx *= rot2(HALFPI * gCurrentMovementProgress);
			}
			break;
		case U:
			if (offset.z == off.z) {
				pa.xy *= rot2(HALFPI * gCurrentMovementProgress);
			}
			break;
		case D:
			if (offset.z == off.y) {
				pa.yx *= rot2(HALFPI * gCurrentMovementProgress);
			}
			break;
		}
		/*
		if (offset.x == off.z) {
			pa.yz *= rot2(gTimeMod);
		}
		if (offset.x == off.y) {
			pa.yz *= rot2(-gTimeMod);
		}
		if (offset.x == off.x) {
			pa.yz *= rot2(-gTimeMod * 0.8);
		}
		if (offset.y == off.z) {
			//pa.xz *= rot2(gTimeMod);
		}
		if (
			!(offset.x == 0. && offset.y == 0.) &&
			!(offset.x == 0. && offset.z == 0.) &&
			!(offset.y == 0. && offset.z == 0.)
		) {
			offset *= sin(gTimeMod) + 2.;
		}
		*/
		offset *= gOffsetStuff;
		vec3 rot = gCubeRot[i];
		if (gCubeCol[1][i] == _x_) {
			cub = centerCube(pa, i, offset, rot);
		} else if (gCubeCol[2][i] == _x_) {
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
		//p.xy *= rot2(gTimeMod/2.);
		//p.yz *= rot2(gTimeMod/3.);
		//p.xz*=rot2(sin(p.z*0.2)*0.2+gTimeMod);
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
	gExclusiveCube = gHitIndex;
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

	// we could not do ambient occlusion and save map() call for more performance
	//float ambientOcc = clamp(map(gHitPosition + normal * .05).x / .05, 0., 1.);
	float ambientOcc = .9; // .9 on purpose, because it's brighter and looks slighty better

	float scat = smoothstep(0., 1., map(gHitPosition - rd * .4).x / .4); // "sub surface scattering"
	shade = mix(specular + shade * (ambientOcc + .2) * (diffuse + scat * .1), shade, fresnel);
	//shade = mix(background, shade, exp(-.002 * result.y * result.y * result.y));
	gExclusiveCube = -1;
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
	gTimeMod = mod(iTime, 10.);
	gCubeCol[0][ 0] = RED; gCubeCol[1][ 0] = YLW; gCubeCol[2][ 0] = BLU;
	gCubeCol[0][ 1] = RED; gCubeCol[1][ 1] = BLU; gCubeCol[2][ 1] = _x_;
	gCubeCol[0][ 2] = RED; gCubeCol[1][ 2] = BLU; gCubeCol[2][ 2] = WHI;
	gCubeCol[0][ 3] = RED; gCubeCol[1][ 3] = YLW; gCubeCol[2][ 3] = _x_;
	gCubeCol[0][ 4] = RED; gCubeCol[1][ 4] = _x_; gCubeCol[2][ 4] = _x_;
	gCubeCol[0][ 5] = RED; gCubeCol[1][ 5] = WHI; gCubeCol[2][ 5] = _x_;
	gCubeCol[0][ 6] = RED; gCubeCol[1][ 6] = GRN; gCubeCol[2][ 6] = YLW;
	gCubeCol[0][ 7] = RED; gCubeCol[1][ 7] = GRN; gCubeCol[2][ 7] = _x_;
	gCubeCol[0][ 8] = RED; gCubeCol[1][ 8] = WHI; gCubeCol[2][ 8] = GRN;
	gCubeCol[0][ 9] = YLW; gCubeCol[1][ 9] = BLU; gCubeCol[2][ 9] = _x_;
	gCubeCol[0][10] = BLU; gCubeCol[1][10] = _x_; gCubeCol[2][10] = _x_;
	gCubeCol[0][11] = BLU; gCubeCol[1][11] = WHI; gCubeCol[2][11] = _x_;
	gCubeCol[0][12] = YLW; gCubeCol[1][12] = _x_; gCubeCol[2][12] = _x_;
	gCubeCol[0][13] = WHI; gCubeCol[1][13] = _x_; gCubeCol[2][13] = _x_;
	gCubeCol[0][14] = GRN; gCubeCol[1][14] = YLW; gCubeCol[2][14] = _x_;
	gCubeCol[0][15] = GRN; gCubeCol[1][15] = _x_; gCubeCol[2][15] = _x_;
	gCubeCol[0][16] = WHI; gCubeCol[1][16] = GRN; gCubeCol[2][16] = _x_;
	gCubeCol[0][17] = YLW; gCubeCol[1][17] = ORG; gCubeCol[2][17] = BLU;
	gCubeCol[0][18] = BLU; gCubeCol[1][18] = ORG; gCubeCol[2][18] = _x_;
	gCubeCol[0][19] = BLU; gCubeCol[1][19] = ORG; gCubeCol[2][19] = WHI;
	gCubeCol[0][20] = YLW; gCubeCol[1][20] = ORG; gCubeCol[2][20] = _x_;
	gCubeCol[0][21] = ORG; gCubeCol[1][21] = _x_; gCubeCol[2][21] = _x_;
	gCubeCol[0][22] = WHI; gCubeCol[1][22] = ORG; gCubeCol[2][22] = _x_;
	gCubeCol[0][23] = GRN; gCubeCol[1][23] = ORG; gCubeCol[2][23] = YLW;
	gCubeCol[0][24] = ORG; gCubeCol[1][24] = GRN; gCubeCol[2][24] = _x_;
	gCubeCol[0][25] = ORG; gCubeCol[1][25] = GRN; gCubeCol[2][25] = WHI;
	gCubeHidden[0] = false;
	gCubeOpacity[0] = 1.;
	gHackFadeStuff = false;
	for (i = 0; i < 26 - 1; i++) {
		int index = gCubeHiddenOrder[i];
		gCubeOpacity[index] = 1.;
		float time = gTimeMod > 9. ? 0. : gTimeMod;
		int whatever = i >= 19 ? 21 : i + 1;
		float until = gNumMovements * MOVEMENT_TIME_SECONDS + float(whatever) * HIDE_TIME_SECONDS;
		if (float(until) < time) {
			gCubeHidden[index] = true;
		} else {
			gCubeHidden[index] = false;
			if (float(until - FADE_TIME_SECONDS) < time) {
				if (i > 21) {
					gHackFadeStuff = true;
				}
				gCubeOpacity[index] = (float(until) - time) / FADE_TIME_SECONDS;
			}
		}
	}

	gCurrentMovement = -1;
	gCurrentMovementProgress = 0.;
	for (i = 0; gTimeMod <= 9. && i < gNumMovements; i++) {
		float until = float(i + 1) * MOVEMENT_TIME_SECONDS;
		if (float(until) < gTimeMod) {
			int tmp;
#define swap(a,b,c,d,e,f,g,h) tmp=gCubeCol[a][b];gCubeCol[a][b]=gCubeCol[c][d];gCubeCol[c][d]=gCubeCol[e][f];gCubeCol[e][f]=gCubeCol[g][h];gCubeCol[g][h]=tmp;
			switch (gMovements[i]) {
			case F:
				swap(0, 2, 1, 0, 1, 17, 2, 19);
				swap(1, 11, 0, 1, 0, 9, 1, 18);
				swap(1, 19, 2, 2, 0, 0, 0, 17);
				swap(1, 2, 2, 0, 2, 17, 0, 19);
				swap(1, 1, 1, 9, 0, 18, 0, 11);
				break;
			case L:
				swap(0, 0, 1, 6, 1, 23, 2, 17);
				swap(2, 0, 0, 6, 0, 23, 1, 17);
				swap(1, 9, 0, 3, 0, 14, 1, 20);
				swap(2, 6, 2, 23, 0, 17, 1, 0);
				swap(1, 3, 1, 14, 0, 20, 0, 9);
				break;
			case R:
				swap(0, 2, 0, 19, 0, 25, 2, 8);
				swap(1, 25, 0, 8, 1, 2, 1, 19);
				swap(0, 5, 0, 11, 1, 22, 1, 16);
				swap(2, 2, 2, 19, 2, 25, 1, 8);
				swap(1, 5, 1, 11, 0, 22, 0, 16);
				break;
			case B:
				swap(0, 8, 2, 25, 1, 23, 2, 6);
				swap(0, 6, 1, 8, 0, 25, 2, 23);
				swap(0, 7, 0, 16, 0, 24, 1, 14);
				swap(2, 8, 1, 25, 0, 23, 1, 6);
				swap(1, 7, 1, 16, 1, 24, 0, 14);
				break;
			case U:
				swap(1, 6, 1, 0, 1, 2, 1, 8);
				swap(2, 8, 2, 6, 2, 0, 2, 2);
				swap(1, 7, 1, 3, 1, 1, 1, 5);
				swap(0, 6, 0, 0, 0, 2, 0, 8);
				swap(0, 7, 0, 3, 0, 1, 0, 5);
				break;
			case D:
				swap(2, 17, 2, 23, 1, 25, 2, 19);
				swap(2, 25, 0, 19, 0, 17, 0, 23);
				swap(0, 18, 0, 20, 1, 24, 0, 22);
				swap(1, 17, 1, 23, 0, 25, 1, 19);
				swap(1, 18, 1, 20, 0, 24, 1, 22);
				break;
			}
		} else {
			gCurrentMovement = gMovements[i];
			gCurrentMovementProgress = 1. - (float(until) - gTimeMod) / MOVEMENT_TIME_SECONDS;
			break;
		}
	}

	vec2 uv=v;uv.y/=1.77;

	vec3 ro = vec3(-80, 80, -70);
	vec3 at = vec3(0, 0, 10);

#if debugmov //noexport
	ro = debug[0].xyz; //noexport
	float vertAngle = debug[1].y/20.; //noexport
	float horzAngle = debug[1].x/20.; //noexport
	if (abs(vertAngle) < .001) { //noexport
		vertAngle = .001; //noexport
	} //noexport
	float xylen = sin(vertAngle); //noexport
	vertAngle = cos(vertAngle); //noexport
	at.x = ro.x + cos(horzAngle) * xylen; //noexport
	at.y = ro.y + sin(horzAngle) * xylen; //noexport
	at.z = ro.z + vertAngle; //noexport
#endif //noexport

	gRounding = ROUNDING;
	gSide = SIDE;
	gUnit = UNIT;
	gOffsetStuff = 1.;
	gShaft = true;
	if (gTimeMod >= 9.) {
		float tt = gTimeMod - 9. / 1.;
		gRounding = mix(3.4, ROUNDING, tt);
		gSide = 12. - gRounding;
		gUnit = gSide * 2. + 2. * gRounding;
		gOffsetStuff = mix(.23, 1., tt);
		ro.z += mix(4., 0., tt);
		at.z += mix(4., 0., tt);
		gShaft = false;
	}

	vec3	cf = normalize(at-ro),
		cl = normalize(cross(cf,vec3(0,0,-1))),
		rd = mat3(cl,normalize(cross(cl,cf)),cf)*normalize(vec3(uv,1)),
		col = vec3(.1) - length(uv) * .05;

	vec4 result = march(ro, rd, 100);

	if (result.x > 0.) { // hit
		vec3 shade = colorHit(result, rd);
		if (gCubeOpacity[gHitIndex] < 1.) {
			float opacity = gCubeOpacity[gHitIndex];
			gCubeHidden[gHitIndex] = true;
			result = march(gHitPosition, rd, 50); // TODO: how many steps?
			vec3 without = result.x > 0. && (!gHackFadeStuff || gHitIndex == 0) ? colorHit(result, rd) : col;
			shade = mix(without, shade, opacity);
		}
		col = shade;
	}

	c = vec4(pow(col, vec3(.4545)), 1.0); // pow for gamma correction because all the cool kids do it
}
