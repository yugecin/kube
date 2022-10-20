// vi: syntax=c
#version 430
#define iTime fpar[0].x
#define debugmov 1 //noexport
#define PI 3.14159265359
#define HALFPI 1.5707963268
#define fight(a,b) if(b.x<a.x)a=b; // union of 2 distance function resuls where x is distance and y is material id
#define BLK 0
#define RED 1
#define BLU 2
#define GRN 3
#define YLW 4
#define WHI 5
#define ORG 6
#define SHA 7
layout (location=0) uniform vec4 fpar[2];
layout (location=2) uniform vec4 debug[2]; //noexport
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

float ROUNDING = 2.5;
float SIDE = 12 - ROUNDING;
float UNIT = SIDE * 2 + 2 * ROUNDING;
float SPACING = UNIT;
vec3 offs = vec3(0., SPACING, -SPACING);
vec4 rot = vec4(0, HALFPI, -HALFPI, PI);

float su(float d1, float d2, float k)
{
	float h = clamp(.5+.5*(d2-d1)/k,0.,1.);
	return mix(d2,d1,h)-k*h*(1.-h);
}
vec3 hitPosition = vec3(0);

vec2 oneSidedCube(vec3 p, int materialTop)
{
	vec2	mc = vec2(length(max(abs(p) - SIDE, 0.)) - ROUNDING, 0),
		tc = vec2(length(max(abs(p + vec3(0., 0., ROUNDING + .02)) - SIDE, 0.)), materialTop);

	return mc.x < tc.x ? mc : tc;
}

vec2 twoSidedCube(vec3 p, int materialTop, int materialFront)
{
	vec2	mc = oneSidedCube(p, materialTop),
		fc = vec2(length(max(abs(p + vec3(0, -ROUNDING - .02, 0.)) - SIDE, 0.)), materialFront);

	return mc.x < fc.x ? mc : fc;
}

// also includes the shaft
vec2 centerCube(vec3 p, int materialTop, vec3 pos, vec3 rot)
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

vec2 middleCube(vec3 p, int materialTop, int materialFront, vec3 pos, vec3 rot)
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

vec2 cornerCube(vec3 p, int materialTop, int materialFront, int materialSide, vec3 pos, vec3 rot)
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
	vec2	aaa = cornerCube(p, RED, YLW, BLU, offs.zyz, rot.yxx),
		aab = middleCube(p, RED, BLU, offs.xyz, rot.xxx),
		aac = cornerCube(p, RED, BLU, WHI, offs.yyz, rot.xxx),
		baa = middleCube(p, RED, YLW, offs.zxz, rot.yxx),
		bab = centerCube(p, RED, offs.xxz, rot.xxx),
		bac = middleCube(p, RED, WHI, offs.yxz, rot.zxx),
		caa = cornerCube(p, RED, GRN, YLW, offs.zzz, rot.wxx),
		cab = middleCube(p, RED, GRN, offs.xzz, rot.wxx),
		cac = cornerCube(p, RED, WHI, GRN, offs.yzz, rot.zxx),
		aba = middleCube(p, YLW, BLU, offs.zyx, rot.xxz),
		abb = centerCube(p, BLU, offs.xyx, rot.xyx),
		abc = middleCube(p, BLU, WHI, offs.yyx, rot.zxz),
		bba = centerCube(p, YLW, offs.zxx, rot.yyx),
		bbc = centerCube(p, WHI, offs.yxx, rot.zyx),
		cba = middleCube(p, GRN, YLW, offs.zzx, rot.yxz),
		cbb = centerCube(p, GRN, offs.xzx, rot.xzx),
		cbc = middleCube(p, WHI, GRN, offs.yzx, rot.wxz),
		aca = cornerCube(p, YLW, ORG, BLU, offs.zyy, rot.yyx),
		acb = middleCube(p, BLU, ORG, offs.xyy, rot.xyx),
		acc = cornerCube(p, BLU, ORG, WHI, offs.yyy, rot.xyx),
		bca = middleCube(p, YLW, ORG, offs.zxy, rot.yyx),
		bcb = centerCube(p, ORG, offs.xxy, rot.xwx),
		bcc = middleCube(p, WHI, ORG, offs.yxy, rot.zyx),
		cca = cornerCube(p, GRN, ORG, YLW, offs.zzy, rot.wyx),
		ccb = middleCube(p, ORG, GRN, offs.xzy, rot.xwx),
		ccc = cornerCube(p, ORG, GRN, WHI, offs.yzy, rot.xwx);

	fight(aaa, aab);
	fight(aaa, aac);
	fight(aaa, baa);
	fight(aaa, bab);
	fight(aaa, bac);
	fight(aaa, caa);
	fight(aaa, cab);
	fight(aaa, cac);
	fight(aaa, aba);
	fight(aaa, abb);
	fight(aaa, abc);
	fight(aaa, bba);
	fight(aaa, bbc);
	fight(aaa, cba);
	fight(aaa, cbb);
	fight(aaa, cbc);
	fight(aaa, aca);
	fight(aaa, acb);
	fight(aaa, acc);
	fight(aaa, bca);
	fight(aaa, bcb);
	fight(aaa, bcc);
	fight(aaa, cca);
	fight(aaa, ccb);
	fight(aaa, ccc);
	return aaa;
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
		if (distance < .03) {
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
	vec2 uv=v;uv.y/=1.77;

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


	vec3	cf = normalize(at-ro),
		cl = normalize(cross(cf,vec3(0,0,-1))),
		rd = mat3(cl,normalize(cross(cl,cf)),cf)*normalize(vec3(uv,1)),
		col = vec3(.1) - length(uv) * .05;

	vec4 result = march(ro, rd, 300);

	if (result.x > 0) { // hit
		vec3 shade = vec3(0);
		float material = result.w;
		if (material == 0) {
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
