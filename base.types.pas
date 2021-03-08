unit base.types;

interface

  uses System.Math;

type


  T_RealArr = packed array of real;



  T_Vect = packed record
    x: double;
    y: double;
    z: double;
    t: double;

  public
    class operator Subtract(const A, B: T_Vect): T_Vect; inline;
    class operator Add(const A, B: T_Vect): T_Vect; inline;
    class operator Divide(const A: T_Vect; const B: real): T_Vect;
    class operator Multiply(const A, B: T_Vect): real; inline;
    class operator Multiply(const A: T_Vect; const b: real): T_Vect; inline;
    class operator Multiply(const A: real; const b: T_Vect): T_Vect; inline;
    class operator Negative(const A: T_Vect): T_Vect;
  end;
  P_Vect = ^T_Vect;

  T_VectArr = array of T_Vect;



  T_Vect4 = array [0..3] of real;

  T_Tens = packed record
    x: T_Vect;
    y: T_Vect;
    z: T_Vect;
  public
    class operator Multiply(const A, B: T_Tens): T_Tens; overload;
    class operator Add(const A: T_Tens; const B: T_Tens): T_Tens; overload;
    class operator Multiply(const A: T_Tens; const B: T_Vect): T_Vect; overload;
    class operator Multiply(const A: T_Vect; const B: T_Tens): T_Vect; overload;
  end;

  P_Tens = ^T_Tens;



  T_Tens3 =  packed record
    xx: double;
    xy: double;
    xz: double;

    yx: double;
    yy: double;
    yz: double;

    zx: double;
    zy: double;
    zz: double;

  public
    function ToTens: T_Tens; overload;
    class operator Explicit(A: T_Tens3): T_Tens; overload;
    class operator Implicit(A: T_Tens): T_Tens3; overload;
  end;



  T_M4 = array [0..3,0..3] of real;
  T_M4Arr = array of T_M4;


  T_TensorArr = array of T_Tens;
  T_TensorArr3 = array of T_Tens3;


  T_SSE = class(TObject)
  private
  public
    class procedure Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j: T_Vect;
        const N: integer); overload; static; register;
    class procedure Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j: T_Vect);
        overload; static; register;
    class procedure Add(var a_j: T_Vect; const C_i: T_Vect; const b_ij: T_Tens);
        overload; static; register;
    class procedure Add(var a_j: T_Vect; const C_i: T_VectArr; const b_ij:
        T_TensorArr; const N: integer); overload; static; register;
    class procedure Add(var a: T_Vect; const C: real; const b: T_Vect; const d:
        real); overload; static; register;
    class procedure Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j:
        T_Vect4); overload; static; register;
    class procedure Add(var a: T_Vect; const b: T_Vect); overload; static;
        register;
    class procedure Add(var a: T_Vect; const b: T_Vect; const C: real); overload;
        static; register;
    class procedure Add(var a_i: T_Vect; const C_j: T_Vect; const b_ij: T_Tens;
        const K: real); overload; static; register;
    class procedure Add(var a: T_Vect; const b : T_Vect; const c, V: real);
        overload; static; register;
    class procedure Add(var a: T_Vect; const C: real; const b: T_Vect); overload;
        static; register;
    class procedure Add(var a_ij: T_M4; const b_ij:
        T_M4; const c: real);         overload;  static;register;
    class procedure Add(var a_ij: T_M4; const DC_i, DC_j:
        T_Vect4);         overload;  static;register;
    class procedure Add(var DC_i: T_Vect4; const a_ij: T_M4;
        const DC_j: T_Vect4); overload; static; register;
    class procedure Add(var a: T_Vect4; const b: T_Vect4); overload;
        static; register;
    class procedure Add(var a: T_Vect4; const b: T_Vect4; const koeff:
        real); overload; static; register;
    class procedure Add(var a: T_Vect4; const b: T_Vect; const koeff: real);
        overload; static; register;
    class procedure Add(var a: T_Tens; const b: T_Tens); overload; static; register;
    class procedure Add(var a: T_Tens; const b: T_Tens; const koeff: real); overload;
        static; register;
    class procedure Invert(var S: T_Tens); overload; static;
    class procedure Add(var a: T_Tens; const b: real); overload; static; register;
    class procedure Add(var B_ab: T_Tens; const _DR_alfa, _DR_beta: T_Vect);
        overload; static; register;
    class procedure Add(var B_ab: T_Tens; const _DR_alfa, _DR_beta: T_Vect; const
        alfa: real); overload; static; register;
    class procedure Add(const a_j: T_VectArr; const C_i: T_VectArr; const b_ij:
        T_TensorArr; const N: integer); overload; static; register;
    class procedure Invert(const S: T_Tens; const i0, i1: integer);
                overload;  static;register;
    class procedure Invert_gauss(var S: T_M4); overload; static; register;
    class procedure Invert_gauss(var S: T_Tens);         overload;  static;register;
    class procedure Invert_gauss(var S: T_M4; const N: integer); overload; static;
        register;
    class function Invert_M3_clang(const M3: T_Tens3): T_Tens3; static; cdecl;
    class procedure Solve(const A: T_Tens; const B: T_Vect; var X: T_Vect);
        overload; static; register;
    class procedure Solve(const A: P_Tens; const B: P_Vect; const X: P_Vect; const
        N: integer); overload; static; register;
    class procedure Solve_old(const A: T_Tens; const B: T_Vect; var X: T_Vect);
        overload; static; register;
    class function DotProd(const V0,V1 : T_Vect): real; overload; static;
    class procedure DotProd(const V0,V1 : T_Vect; var Res: double); overload; static;
    class procedure DotProd(const V0,V1 : T_VectArr;const Res: T_RealArr; const i0,i1:
        integer); overload; static;
    class procedure DotProd(const V0,V1 : T_VectArr; var Res: real; const i0,i1:
        integer); overload; static;
  end;


procedure Add(var a: real; const b : real); overload; register; inline;

procedure Add(var a: real; const b, c : T_Vect); overload; register; inline;

procedure Add(var a: T_M4; const b: real); overload; inline; register;

procedure Add(var a_ij: T_M4; const b_ij: T_M4);
    overload; inline; register;

procedure Add(var a_ij: T_M4; const b_ij: T_M4;
    const c: real); overload; inline; register;

procedure Add(var a_ij: T_M4; const  DC_i, DC_j: T_Vect4);
    overload; inline; register;

procedure Add(var DC_i: T_Vect4; const a_ij: T_M4; const
    DC_j: T_Vect4); overload; inline; register;

procedure Add(var a: T_Vect4; const b: T_Vect4); overload; inline; register;

procedure Add(var a: T_Vect4; const b: T_Vect4; const koeff: real);
    overload; inline; register;

procedure Add(var a: T_Tens; const b: T_Tens); overload; inline; register;

procedure Add(var a: T_Tens; const b: T_Tens; const koeff: real); overload; inline; register;

procedure Add(var a: T_Tens; const b: real); overload; inline; register;

procedure Add(var a: T_Vect; const b: T_Vect; const C: real); overload;
    register; inline;

procedure Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j: T_Vect4);
    overload; inline; register;

procedure Add(var a: T_Vect; const b: T_Vect);overload; register; inline;

procedure Add(var a: T_Vect; const b : T_Vect; const c, V: real); overload;
    register; inline;

procedure Add(var a_i: T_Vect; const C_j: T_Vect; const b_ij: T_Tens; const K:
    real); overload; register; inline;

function CubeRoot(A: real): real;

function Det(const mm: T_Tens; const a : real): real; inline;

function Divide(const r: T_Vect; const h: T_Tens): T_Vect; overload; inline;

function EQ(const a,b : T_Vect; const eps : real = 1.0e-11): Boolean; inline;

function Matrix3(const DiagonalVal : real): T_Tens; overload; inline;

function Matrix3(const t00,t01,t02, t10,t11,t12, t20,t21,t22 : real): T_Tens;
    overload; inline;

function Matrix4(const DiagonalVal : real): T_M4; overload;
    inline;

function Matrix4(const t00,t01,t02,t03, t10,t11,t12,t13, t20,t21,t22,t23, t30,
    t31,t32,t33 : real): T_M4; overload; inline;

function mult(const M: array of T_M4): T_M4; overload;

function mult(const V0, V1: real): real; overload; inline;

function mult(const Mik, Mkj : T_M4): T_M4;
    overload; inline;

function mult(const Mik: T_M4; const Pk: T_Vect4): T_Vect4; overload; inline;

function mult(const P: T_Vect4; const a: real): T_Vect4; overload;
    inline;

function mult(const Pk: T_Vect4; const Mki: T_M4): T_Vect4; overload; inline;

function mult(const M1: T_M4; const S : real): T_M4; overload; inline;

function mult(const Mik : T_Tens; const c: real): T_Tens; overload; inline;

function mult(const Mij: T_Tens; const Vj: T_Vect): T_Vect; overload; inline;

function mult_transponed(const Mik: T_Tens; const Mjk : T_Tens): T_Tens; overload;
    inline;





function mult(const a : T_vect; const Lmd : real): T_Vect; overload; register;
    inline;

function mult(const Pk: T_Vect; const Mki : T_M4): T_Vect; overload; inline;

function mult(const Vk : T_Vect;const Mki : T_Tens): T_Vect; overload; inline;



function mult(const V0,V1 : T_Vect): T_Vect; overload; inline;

function Norma(const DC : T_Vect4): real; overload; inline;

function Norma(const S : T_M4): real; overload; inline;

function Norma(const S : T_M4; const Diag : real): real; overload; inline;

function Norma(const S1,S2 : T_M4): real; overload; inline;

function Norma(const S : T_Tens): real; overload; inline;

function Norma(const S1,S2 : T_Tens): real; overload; inline;

function Norma(const V : T_Vect): real; overload; inline;

function Norma(const V1,V2 : T_Vect): real; overload; inline;

function razn(const a,b : real): real; overload; inline; register;

function razn(const a,b : T_Tens): T_Tens; overload; inline; register;

function razn(const a,b : T_Vect): T_Vect; overload; inline; register;

procedure SetUnite(out a : T_M4); overload; inline;

procedure set_bit(var b: byte; const bit: byte; const st: boolean); inline;

function summ(const a: array of T_Vect): T_Vect; overload;

function summ(const a: T_Tens; const b : real): T_Tens; overload; inline;

function summ(const a,b : T_Tens): T_Tens; overload; inline;

function summ(const a: T_Vect; const b : real): T_Vect; overload; inline;

function summ(const a,b : T_Vect): T_Vect; overload; inline;

function summ_t(const a: array of T_Tens): T_Tens;

function svertka(const M0, M1 : T_Tens): real; overload; inline;

function tensMax(const a,b : T_Tens): T_Tens;

function Tensor(const DiagonalVal : real): T_Tens; overload; inline;

function Tensor(const t00,t01,t02, t10,t11,t12, t20,t21,t22 : real): T_Tens;
    overload; inline;

function Tensor(const F0,F1,F2 : T_Vect): T_Tens; overload; inline;

procedure Transpone(var S: T_M4); overload; inline;

procedure Transpone(var S: T_Tens); overload; inline;

function Vect(const X,Y,Z : real): T_Vect; overload; inline;

function Vect(const DC : T_Vect4): T_Vect; overload; inline;

function VectCube(const a : T_Vect): real; overload;

function VectCube(const a,b : T_Vect): real; overload;

function Vector4(const v0,v1,v2,v3: real): T_Vect4; overload; inline;

function Volume(const S: T_Tens): real; overload; inline;

function Volume(const S0,S1,S2: T_Vect): real; overload; inline;

function Volume(const center, S0,S1,S2: T_Vect): real; overload; inline;

function VSqr(const M4: T_M4): real; overload; inline;

function VSqr(const a,b : T_Tens): real; overload; inline;

function VSqr(const a: T_Vect): real; overload; register; inline;

function VSqr(const a,b: T_Vect): real; overload; register; inline;

function DotProd(const V0,V1 : T_Vect4): real; overload; inline;

function DotProd(const V0,V1 : T_Vect): real; overload; register; inline;

procedure Zero(out a : integer); overload; register; inline;

procedure Zero(out a : T_M4); overload; register; inline;

procedure Zero(out a : T_Vect4); overload; register; inline;

procedure Zero(out a : T_RealArr); overload; register; inline;

procedure Zero(out a : T_Tens); overload; register; inline;

procedure Zero(out a : T_TensorArr); overload; inline;

procedure Zero(out a : T_Vect); overload; register; inline;

procedure Zero(out a : T_VectArr); overload; register; inline;

procedure swap(var a,b: real); register; inline;

procedure Zero(out a : real); overload; register; inline;

function Norma(const V : T_RealArr; const cnt: integer = -1): real; overload;

procedure MinVect(var a: T_Vect; const b: T_Vect);

procedure MaxVect(var a: T_Vect; const b: T_Vect);

const m_pi : extended = 3.1415926535897932384626433832795;

const Tens0 : T_Tens =  ( X:( X:0.0; Y:0.0; Z:0.0; t:0.0; ) ;
                          Y:( X:0.0; Y:0.0; Z:0.0; t:0.0; ) ;
                          Z:( X:0.0; Y:0.0; Z:0.0; t:0.0; ) );
const Vect0 : T_Vect = (X:0.0;Y:0.0;Z:0.0;t:0.0);
const S1: T_Tens =      ( X:( X:1.0; Y:0.0; Z:0.0; t:0.0; ) ;
                          Y:( X:0.0; Y:1.0; Z:0.0; t:0.0; ) ;
                          Z:( X:0.0; Y:0.0; Z:1.0; t:0.0; ) );

const M4_Unite: T_M4 = ( ( 1.0, 0.0, 0.0, 0.0 ),
                         ( 0.0, 1.0, 0.0, 0.0 ),
                         ( 0.0, 0.0, 1.0, 0.0 ),
                         ( 0.0, 0.0, 0.0, 1.0 ));

const M4_ZERO: T_M4 = (   (  0.0,  0.0,  0.0,  0.0),
                          (  0.0,  0.0,  0.0,  0.0),
                          (  0.0,  0.0,  0.0,  0.0),
                          (  0.0,  0.0,  0.0,  0.0));

  function ALIGN(const N: integer): integer; overload; inline;
  function ALIGN(const N, Step: integer): integer; overload; inline;

function Norma(const V1,V2 : T_RealArr; const cnt: integer = -1): real;
    overload;

procedure Add(var a: T_Vect; const C: real; const b: T_Vect; const d: real);
    overload; register; inline;

procedure Add(var B_ab: T_Tens; const _DR_alfa, _DR_beta: T_Vect; const alfa:
    real); overload; inline;

function mult(const Mij: T_Tens; const Vj, Vi: T_Vect): real; overload; inline;

procedure Add(var a: T_Vect; const C: real; const b: T_Vect); overload; inline;  register;

procedure SetUnite(out a: T_Tens); overload; inline;

procedure Add(var a: real; const b, c : T_Vect; const V: real); overload;
    inline;    register;

procedure Add(var B_ab: T_Tens; const _DR_alfa, _DR_beta: T_Vect); overload;
    inline;

procedure Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j: T_Vect);
    overload; register; inline;

procedure Invert(var S: T_Tens3); overload;

procedure Invert_delphi_assembler_cleared(var S: T_Tens);

procedure Invert_gauss(var S: T_Tens);

procedure Invert(var a : T_M4); overload; inline;

procedure Invert(const a : T_M4 ; out E : T_M4); overload; inline;

procedure Gauss(const Matr: T_Tens; const FreePart: T_Vect; var Res: T_Vect);

function Vector3(const v0,v1,v2: real): T_Vect; overload; inline;

function VSqr(const a : T_Tens): real; overload; inline;

procedure Invert_m4(var a: T_M4); inline;

procedure Invert2(var S: T_Tens);

procedure Invert(var S: T_Tens); overload;

function mult(const Mik, Mkj : T_Tens): T_Tens; overload; inline;

implementation

const FloatSize = 8;
const VectSize = FloatSize*4;



  function ALIGN(const N: integer): integer; overload; inline;
  begin
    result := N + 8 - (N mod 8);
  end;
  function ALIGN(const N, Step: integer): integer; overload; inline;
  begin
    result := N + Step - (N mod Step);
  end;


procedure Add(var a: real; const b : real);
begin
  a := a+b;
end;

procedure Add(var a: real; const b, c : T_Vect);
begin
  a := a + DotProd( b , c );
end;

procedure Add(var a: T_M4; const b: real);
begin
  a[0,0] := a[0,0] + b;
  a[1,1] := a[1,1] + b;
  a[2,2] := a[2,2] + b;
  a[3,3] := a[3,3] + b;
end;

procedure Add(var a_ij: T_M4; const b_ij: T_M4);
begin
  a_ij[ 0, 0] := a_ij[ 0, 0] + b_ij[ 0, 0];
  a_ij[ 0, 1] := a_ij[ 0, 1] + b_ij[ 0, 1];
  a_ij[ 0, 2] := a_ij[ 0, 2] + b_ij[ 0, 2];
  a_ij[ 0, 3] := a_ij[ 0, 3] + b_ij[ 0, 3];

  a_ij[ 1, 0] := a_ij[ 1, 0] + b_ij[ 1, 0];
  a_ij[ 1, 1] := a_ij[ 1, 1] + b_ij[ 1, 1];
  a_ij[ 1, 2] := a_ij[ 1, 2] + b_ij[ 1, 2];
  a_ij[ 1, 3] := a_ij[ 1, 3] + b_ij[ 1, 3];

  a_ij[ 2, 0] := a_ij[ 2, 0] + b_ij[ 2, 0];
  a_ij[ 2, 1] := a_ij[ 2, 1] + b_ij[ 2, 1];
  a_ij[ 2, 2] := a_ij[ 2, 2] + b_ij[ 2, 2];
  a_ij[ 2, 3] := a_ij[ 2, 3] + b_ij[ 2, 3];

  a_ij[ 3, 0] := a_ij[ 3, 0] + b_ij[ 3, 0];
  a_ij[ 3, 1] := a_ij[ 3, 1] + b_ij[ 3, 1];
  a_ij[ 3, 2] := a_ij[ 3, 2] + b_ij[ 3, 2];
  a_ij[ 3, 3] := a_ij[ 3, 3] + b_ij[ 3, 3];
end;

procedure Add(var a_ij: T_M4; const b_ij: T_M4;
    const c: real);
begin
  a_ij[ 0, 0] := a_ij[ 0, 0] + b_ij[ 0, 0]*c;
  a_ij[ 0, 1] := a_ij[ 0, 1] + b_ij[ 0, 1]*c;
  a_ij[ 0, 2] := a_ij[ 0, 2] + b_ij[ 0, 2]*c;
  a_ij[ 0, 3] := a_ij[ 0, 3] + b_ij[ 0, 3]*c;

  a_ij[ 1, 0] := a_ij[ 1, 0] + b_ij[ 1, 0]*c;
  a_ij[ 1, 1] := a_ij[ 1, 1] + b_ij[ 1, 1]*c;
  a_ij[ 1, 2] := a_ij[ 1, 2] + b_ij[ 1, 2]*c;
  a_ij[ 1, 3] := a_ij[ 1, 3] + b_ij[ 1, 3]*c;

  a_ij[ 2, 0] := a_ij[ 2, 0] + b_ij[ 2, 0]*c;
  a_ij[ 2, 1] := a_ij[ 2, 1] + b_ij[ 2, 1]*c;
  a_ij[ 2, 2] := a_ij[ 2, 2] + b_ij[ 2, 2]*c;
  a_ij[ 2, 3] := a_ij[ 2, 3] + b_ij[ 2, 3]*c;

  a_ij[ 3, 0] := a_ij[ 3, 0] + b_ij[ 3, 0]*c;
  a_ij[ 3, 1] := a_ij[ 3, 1] + b_ij[ 3, 1]*c;
  a_ij[ 3, 2] := a_ij[ 3, 2] + b_ij[ 3, 2]*c;
  a_ij[ 3, 3] := a_ij[ 3, 3] + b_ij[ 3, 3]*c;
end;

procedure Add(var a_ij: T_M4; const  DC_i, DC_j: T_Vect4);
begin
  a_ij[ 0, 0] := a_ij[ 0, 0] + DC_j[ 0]*DC_i[ 0];
  a_ij[ 0, 1] := a_ij[ 0, 1] + DC_j[ 1]*DC_i[ 0];
  a_ij[ 0, 2] := a_ij[ 0, 2] + DC_j[ 2]*DC_i[ 0];
  a_ij[ 0, 3] := a_ij[ 0, 3] + DC_j[ 3]*DC_i[ 0];

  a_ij[ 1, 0] := a_ij[ 1, 0] + DC_j[ 0]*DC_i[ 1];
  a_ij[ 1, 1] := a_ij[ 1, 1] + DC_j[ 1]*DC_i[ 1];
  a_ij[ 1, 2] := a_ij[ 1, 2] + DC_j[ 2]*DC_i[ 1];
  a_ij[ 1, 3] := a_ij[ 1, 3] + DC_j[ 3]*DC_i[ 1];

  a_ij[ 2, 0] := a_ij[ 2, 0] + DC_j[ 0]*DC_i[ 2];
  a_ij[ 2, 1] := a_ij[ 2, 1] + DC_j[ 1]*DC_i[ 2];
  a_ij[ 2, 2] := a_ij[ 2, 2] + DC_j[ 2]*DC_i[ 2];
  a_ij[ 2, 3] := a_ij[ 2, 3] + DC_j[ 3]*DC_i[ 2];

  a_ij[ 3, 0] := a_ij[ 3, 0] + DC_j[ 0]*DC_i[ 3];
  a_ij[ 3, 1] := a_ij[ 3, 1] + DC_j[ 1]*DC_i[ 3];
  a_ij[ 3, 2] := a_ij[ 3, 2] + DC_j[ 2]*DC_i[ 3];
  a_ij[ 3, 3] := a_ij[ 3, 3] + DC_j[ 3]*DC_i[ 3];

end;

procedure Add(var DC_i: T_Vect4; const a_ij: T_M4; const
    DC_j: T_Vect4);
begin
  DC_i[ 0] := DC_i[ 0] + a_ij[ 0, 0] *DC_j[ 0]+a_ij[ 0, 1] *DC_j[ 1]+a_ij[ 0, 2] *DC_j[ 2]+a_ij[ 0, 3] *DC_j[ 3];
  DC_i[ 1] := DC_i[ 1] + a_ij[ 1, 0] *DC_j[ 0]+a_ij[ 1, 1] *DC_j[ 1]+a_ij[ 1, 2] *DC_j[ 2]+a_ij[ 1, 3] *DC_j[ 3];
  DC_i[ 2] := DC_i[ 2] + a_ij[ 2, 0] *DC_j[ 0]+a_ij[ 2, 1] *DC_j[ 1]+a_ij[ 2, 2] *DC_j[ 2]+a_ij[ 2, 3] *DC_j[ 3];
  DC_i[ 3] := DC_i[ 3] + a_ij[ 3, 0] *DC_j[ 0]+a_ij[ 3, 1] *DC_j[ 1]+a_ij[ 3, 2] *DC_j[ 2]+a_ij[ 3, 3] *DC_j[ 3];
end;

procedure Add(var a: T_Vect4; const b: T_Vect4);
begin
  a[0] := a[0]+b[0];
  a[1] := a[1]+b[1];
  a[2] := a[2]+b[2];
  a[3] := a[3]+b[3];
end;

procedure Add(var a: T_Vect4; const b: T_Vect4; const koeff: real);
begin
  a[ 0] := a[ 0]+b[ 0]*koeff;
  a[ 1] := a[ 1]+b[ 1]*koeff;
  a[ 2] := a[ 2]+b[ 2]*koeff;
  a[ 3] := a[ 3]+b[ 3]*koeff;
end;

procedure Add(var a: T_Tens; const b: T_Tens);
begin
  a.x.x := a.x.x + b.x.x;
  a.y.x := a.y.x + b.y.x;
  a.z.x := a.z.x + b.z.x;

  a.x.y := a.x.y + b.x.y;
  a.y.y := a.y.y + b.y.y;
  a.z.y := a.z.y + b.z.y;

  a.x.z := a.x.z + b.x.z;
  a.y.z := a.y.z + b.y.z;
  a.z.z := a.z.z + b.z.z;
end;

procedure Add(var a: T_Tens; const b: T_Tens; const koeff: real);
begin
  a.x.x := a.x.x + b.x.x*koeff;
  a.x.y := a.x.y + b.x.y*koeff;
  a.x.z := a.x.z + b.x.z*koeff;

  a.y.x := a.y.x + b.y.x*koeff;
  a.y.y := a.y.y + b.y.y*koeff;
  a.y.z := a.y.z + b.y.z*koeff;

  a.z.x := a.z.x + b.z.x*koeff;
  a.z.y := a.z.y + b.z.y*koeff;
  a.z.z := a.z.z + b.z.z*koeff;
end;

procedure Add(var a: T_Tens; const b: real);
begin
  a.x.x := a.x.x + b;
  a.y.y := a.y.y + b;
  a.z.z := a.z.z + b;
end;

procedure Add(var a: T_Vect; const b: T_Vect; const C: real);
begin
  a.x := a.x+b.x*C;
  a.y := a.y+b.y*C;
  a.z := a.z+b.z*C;
end;

procedure Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j: T_Vect4);
begin
  a_i.x := a_i.x+ b_ij.x.x*C_j[0]+b_ij.x.y*C_j[1]+b_ij.x.z*C_j[2] ;
  a_i.y := a_i.y+ b_ij.y.x*C_j[0]+b_ij.y.y*C_j[1]+b_ij.y.z*C_j[2] ;
  a_i.z := a_i.z+ b_ij.z.x*C_j[0]+b_ij.z.y*C_j[1]+b_ij.z.z*C_j[2] ;
end;

procedure Add(var a: T_Vect; const b: T_Vect);
begin
  a.x := a.x+b.x;
  a.y := a.y+b.y;
  a.z := a.z+b.z;
end;

procedure Add(var a: T_Vect; const b : T_Vect; const c, V: real);
begin
  a.x := a.x + b.x*c*V;
  a.y := a.y + b.y*c*V;
  a.z := a.z + b.z*c*V;
end;

procedure Add(var a_i: T_Vect; const C_j: T_Vect; const b_ij: T_Tens; const K:
    real);
begin
  a_i.x := a_i.x+K*(b_ij.x.x*C_j.x+b_ij.x.y*C_j.y+b_ij.x.z*C_j.z) ;
  a_i.y := a_i.y+K*(b_ij.y.x*C_j.x+b_ij.y.y*C_j.y+b_ij.y.z*C_j.z) ;
  a_i.z := a_i.z+K*(b_ij.z.x*C_j.x+b_ij.z.y*C_j.y+b_ij.z.z*C_j.z) ;
end;

function CubeRoot(A: real): real;
var x: real;
begin

  if A < 0  then begin
    result := - CubeRoot(-A);
    exit;
  end;

  if A =0.0 then begin
    result := 0.0;
    exit;
  end;

  if A > 4096 then begin
    result := 16*CubeRoot(A/4096);
    exit;

  end;

  if A > 1.0 then
  begin
    x := sqrt(A);

    x := 0.333333333333333*( 2*x  + A / sqr(x) );
    x := 0.333333333333333*( 2*x  + A / sqr(x) );
    x := 0.333333333333333*( 2*x  + A / sqr(x) );
    x := 0.333333333333333*( 2*x  + A / sqr(x) );

    x := 0.333333333333333*( 2*x  + A / sqr(x) );
    x := 0.333333333333333*( 2*x  + A / sqr(x) );
    x := 0.333333333333333*( 2*x  + A / sqr(x) );
    x := 0.333333333333333*( 2*x  + A / sqr(x) );


    result := x;
    exit;
  end;

  result := 1.0/CubeRoot( 1.0/ A);

end;

function Det(const mm: T_Tens; const a : real): real;
begin
  result :=  ( mm.x.x-a )*( ( mm.y.y-a )*( mm.z.z-a )  -  ( mm.y.z   )*mm.z.y )
            -( mm.y.x   )*( ( mm.x.y   )*( mm.z.z-a )  -  ( mm.z.y   )*mm.x.z )
            +( mm.z.x   )*( ( mm.x.y   )*( mm.y.z   )  -  ( mm.y.y-a )*mm.x.z );
end;

function Divide(const r: T_Vect; const h: T_Tens): T_Vect;
var  D : real;
begin
  D:=    h.x.x*( h.y.y*h.z.z - h.y.z*h.z.y )
       + h.x.y*( h.z.x*h.y.z - h.y.x*h.z.z )
       + h.x.z*( h.y.x*h.z.y - h.y.y*h.z.x ) ;

  if Abs(D)> 1.0e-16 then
    D := 1.0 / D
  else begin
    result := r;
    exit;
  end;


  Result.x :=(    (  h.y.y*h.z.z - h.y.z*h.z.y  ) * r.x
                 -(  h.x.y*h.z.z - h.x.z*h.z.y  ) * r.y
                 +(  h.x.y*h.y.z - h.x.z*h.y.y  ) * r.z  ) * D;

  Result.y :=(  -(  h.y.x*h.z.z - h.z.x*h.y.z  ) * r.x
                 +(  h.x.x*h.z.z - h.x.z*h.z.x  ) * r.y
                 -(  h.x.x*h.y.z - h.y.x*h.x.z  ) * r.z  ) * D;

  Result.z :=(   (  h.y.x*h.z.y - h.y.y*h.z.x  ) * r.x
                 -(  h.x.x*h.z.y - h.x.y*h.z.x  ) * r.y
                 +(  h.x.x*h.y.y - h.x.y*h.y.x  ) * r.z  ) * D;

end;

function EQ(const a,b : T_Vect; const eps : real = 1.0e-11): Boolean;
begin
  Result := ( abs( a.x-b.x ) < eps ) and
            ( abs( a.y-b.y ) < eps ) and
            ( abs( a.z-b.z ) < eps );
end;

procedure Gauss(const Matr: T_Tens; const FreePart: T_Vect; var Res: T_Vect);
  var M: array [0..2, 0..3] of real; // Extended matrix
      t: real;

begin
  M[0, 0] := Matr.x.x;  M[0, 1] := Matr.x.y;  M[0, 2] := Matr.x.z;
  M[1, 0] := Matr.y.x;  M[1, 1] := Matr.y.y;  M[1, 2] := Matr.y.z;
  M[2, 0] := Matr.z.x;  M[2, 1] := Matr.z.y;  M[2, 2] := Matr.z.z;

  M[0, 3] := FreePart.x;
  M[1, 3] := FreePart.y;
  M[2, 3] := FreePart.z;

  // коэффициенты матрицы получились удвоенные (сумма H, вместо полусуммы)
  // значит, итоговое значение dR будет располовиненным, то есть,
  // вычисленное значение u надо будет умножить на 2.0, а детерминант- на 1/8
//  negativeDet := false;
  if (abs(M[2, 0])>abs(M[0, 0])) and (abs(M[2, 0]) > abs(M[1, 0])) then
  begin
    t := M[2, 0]; M[2, 0] := M[0, 0]; M[0, 0]:= t;
    t := M[2, 1]; M[2, 1] := M[0, 1]; M[0, 1]:= t;
    t := M[2, 2]; M[2, 2] := M[0, 2]; M[0, 2]:= t;
    t := M[2, 3]; M[2, 3] := M[0, 3]; M[0, 3]:= t;
  end;

  if abs(M[1, 0]) > abs(M[0, 0]) then
  begin
    t := M[1, 0]; M[1, 0] := M[0, 0]; M[0, 0]:= t;
    t := M[1, 1]; M[1, 1] := M[0, 1]; M[0, 1]:= t;
    t := M[1, 2]; M[1, 2] := M[0, 2]; M[0, 2]:= t;
    t := M[1, 3]; M[1, 3] := M[0, 3]; M[0, 3]:= t;
  end;

  t := M[1, 0]/M[0, 0];
  M[1, 1] := M[1, 1] - t*M[0, 1];
  M[1, 2] := M[1, 2] - t*M[0, 2];
  M[1, 3] := M[1, 3] - t*M[0, 3];

  t := M[2, 0]/M[0, 0];
  M[2, 1] := M[2, 1] - t*M[0, 1];
  M[2, 2] := M[2, 2] - t*M[0, 2];
  M[2, 3] := M[2, 3] - t*M[0, 3];

  if Abs(M[2, 1])>Abs(M[1, 1]) then
  begin
//    t := M[2, 0]; M[2, 0] := M[1, 0]; M[1, 0]:= t;          // Zeroth columnt may not be moved, It fulfilled by zeros
    t := M[2, 1]; M[2, 1] := M[1, 1]; M[1, 1]:= t;
    t := M[2, 2]; M[2, 2] := M[1, 2]; M[1, 2]:= t;
    t := M[2, 3]; M[2, 3] := M[1, 3]; M[1, 3]:= t;
  end;

  t := M[2, 1]/M[1, 1];
  M[2, 2] := M[2, 2]-t*M[1, 2];
  M[2, 3] := M[2, 3]-t*M[1, 3];

  M[2, 3] := M[2, 3]/M[2, 2];
  {x2}

  M[1, 3] := ( M[1, 3]-M[1, 2]*M[2, 3] )/M[1, 1];
  {x1}        {F1}           {x2}

  M[0, 3] := (M[0, 3]- M[0, 2]*M[2, 3] - M[0, 1]* M[1, 3])/M[0, 0];
  {x0}         {F0}          {x2}             {x1}



  Res.x := M[0, 3];
  Res.y := M[1, 3];
  Res.z := M[2, 3];


end;

function Matrix3(const DiagonalVal : real): T_Tens;
begin
  result .x.x := DiagonalVal;
  result .x.y := 0.0;
  result .x.z := 0.0;

  result .y.x := 0.0;
  result .y.y := DiagonalVal;
  result .y.z := 0.0;

  result .z.x := 0.0;
  result .z.y := 0.0;
  result .z.z := DiagonalVal;


end;

function Matrix3(const t00,t01,t02, t10,t11,t12, t20,t21,t22 : real): T_Tens;
begin
  result .x.x := t00;
  result .x.y := t01;
  result .x.z := t02;

  result .y.x := t10;
  result .y.y := t11;
  result .y.z := t12;

  result .z.x := t20;
  result .z.y := t21;
  result .z.z := t22;


end;

function Matrix4(const DiagonalVal : real): T_M4;
begin
  ZERO(Result);

  result [ 0, 0] := DiagonalVal;
  result [ 1, 1] := DiagonalVal;
  result [ 2, 2] := DiagonalVal;
  result [ 3, 3] := DiagonalVal;
end;

function Matrix4(const t00,t01,t02,t03, t10,t11,t12,t13, t20,t21,t22,t23, t30,
    t31,t32,t33 : real): T_M4;
begin
  T_M4(result)[0,0] := t00;
  T_M4(result)[0,1] := t01;
  T_M4(result)[0,2] := t02;
  T_M4(result)[0,3] := t03;

  T_M4(result)[1,0] := t10;
  T_M4(result)[1,1] := t11;
  T_M4(result)[1,2] := t12;
  T_M4(result)[1,3] := t13;


  T_M4(result)[2,0] := t20;
  T_M4(result)[2,1] := t21;
  T_M4(result)[2,2] := t22;
  T_M4(result)[2,3] := t23;


  T_M4(result)[3,0] := t30;
  T_M4(result)[3,1] := t31;
  T_M4(result)[3,2] := t32;
  T_M4(result)[3,3] := t33;


end;

function Max(const a,b : real): real; inline;
begin
  if a>b then Result := a     else Result := b;
end;

function Min(const a,b : real): real; inline;
begin
  if a>b then Result := b     else Result := a;
end;

function mult(const M: array of T_M4): T_M4;
var i : integer;
begin
  result := M[0];
  for i := 1 to high(M) do
    result:= Mult( result, M[i] );

end;

function mult(const V0, V1: real): real;
begin
  Result := V0*V1;
end;

function mult(const M1: T_M4; const S : real):
    T_M4;
var i,j : integer;
begin
  for i := 0 to 3 do
  for j := 0 to 3 do
    result[i,j]:= M1[i,j]*S;
end;

function mult(const Mik, Mkj : T_M4): T_M4;
begin
    result[ 0, 0]:=  Mik[ 0, 0]*Mkj[ 0, 0]+ Mik[ 0, 1]*Mkj[ 1, 0]+ Mik[ 0, 2]*Mkj[ 2, 0]+ Mik[ 0, 3]*Mkj[ 3, 0];
    result[ 0, 1]:=  Mik[ 0, 0]*Mkj[ 0, 1]+ Mik[ 0, 1]*Mkj[ 1, 1]+ Mik[ 0, 2]*Mkj[ 2, 1]+ Mik[ 0, 3]*Mkj[ 3, 1];
    result[ 0, 2]:=  Mik[ 0, 0]*Mkj[ 0, 2]+ Mik[ 0, 1]*Mkj[ 1, 2]+ Mik[ 0, 2]*Mkj[ 2, 2]+ Mik[ 0, 3]*Mkj[ 3, 2];
    result[ 0, 3]:=  Mik[ 0, 0]*Mkj[ 0, 3]+ Mik[ 0, 1]*Mkj[ 1, 3]+ Mik[ 0, 2]*Mkj[ 2, 3]+ Mik[ 0, 3]*Mkj[ 3, 3];

    result[ 1, 0]:=  Mik[ 1, 0]*Mkj[ 0, 0]+ Mik[ 1, 1]*Mkj[ 1, 0]+ Mik[ 1, 2]*Mkj[ 2, 0]+ Mik[ 1, 3]*Mkj[ 3, 0];
    result[ 1, 1]:=  Mik[ 1, 0]*Mkj[ 0, 1]+ Mik[ 1, 1]*Mkj[ 1, 1]+ Mik[ 1, 2]*Mkj[ 2, 1]+ Mik[ 1, 3]*Mkj[ 3, 1];
    result[ 1, 2]:=  Mik[ 1, 0]*Mkj[ 0, 2]+ Mik[ 1, 1]*Mkj[ 1, 2]+ Mik[ 1, 2]*Mkj[ 2, 2]+ Mik[ 1, 3]*Mkj[ 3, 2];
    result[ 1, 3]:=  Mik[ 1, 0]*Mkj[ 0, 3]+ Mik[ 1, 1]*Mkj[ 1, 3]+ Mik[ 1, 2]*Mkj[ 2, 3]+ Mik[ 1, 3]*Mkj[ 3, 3];

    result[ 2, 0]:=  Mik[ 2, 0]*Mkj[ 0, 0]+ Mik[ 2, 1]*Mkj[ 1, 0]+ Mik[ 2, 2]*Mkj[ 2, 0]+ Mik[ 2, 3]*Mkj[ 3, 0];
    result[ 2, 1]:=  Mik[ 2, 0]*Mkj[ 0, 1]+ Mik[ 2, 1]*Mkj[ 1, 1]+ Mik[ 2, 2]*Mkj[ 2, 1]+ Mik[ 2, 3]*Mkj[ 3, 1];
    result[ 2, 2]:=  Mik[ 2, 0]*Mkj[ 0, 2]+ Mik[ 2, 1]*Mkj[ 1, 2]+ Mik[ 2, 2]*Mkj[ 2, 2]+ Mik[ 2, 3]*Mkj[ 3, 2];
    result[ 2, 3]:=  Mik[ 2, 0]*Mkj[ 0, 3]+ Mik[ 2, 1]*Mkj[ 1, 3]+ Mik[ 2, 2]*Mkj[ 2, 3]+ Mik[ 2, 3]*Mkj[ 3, 3];

    result[ 3, 0]:=  Mik[ 3, 0]*Mkj[ 0, 0]+ Mik[ 3, 1]*Mkj[ 1, 0]+ Mik[ 3, 2]*Mkj[ 2, 0]+ Mik[ 3, 3]*Mkj[ 3, 0];
    result[ 3, 1]:=  Mik[ 3, 0]*Mkj[ 0, 1]+ Mik[ 3, 1]*Mkj[ 1, 1]+ Mik[ 3, 2]*Mkj[ 2, 1]+ Mik[ 3, 3]*Mkj[ 3, 1];
    result[ 3, 2]:=  Mik[ 3, 0]*Mkj[ 0, 2]+ Mik[ 3, 1]*Mkj[ 1, 2]+ Mik[ 3, 2]*Mkj[ 2, 2]+ Mik[ 3, 3]*Mkj[ 3, 2];
    result[ 3, 3]:=  Mik[ 3, 0]*Mkj[ 0, 3]+ Mik[ 3, 1]*Mkj[ 1, 3]+ Mik[ 3, 2]*Mkj[ 2, 3]+ Mik[ 3, 3]*Mkj[ 3, 3];

end;

function mult(const P: T_Vect4; const a: real): T_Vect4;
begin
  result[ 0]:= P[ 0]*a;
  result[ 1]:= P[ 1]*a;
  result[ 2]:= P[ 2]*a;
  result[ 3]:= P[ 3]*a;
end;

function mult(const Mik: T_M4; const Pk: T_Vect4): T_Vect4;
begin
  result[0]:= Pk[0]*Mik[0,0]+ Pk[1]*Mik[0,1]+ Pk[2]*Mik[0,2]+ Pk[3]*Mik[0,3];
  result[1]:= Pk[0]*Mik[1,0]+ Pk[1]*Mik[1,1]+ Pk[2]*Mik[1,2]+ Pk[3]*Mik[1,3];
  result[2]:= Pk[0]*Mik[2,0]+ Pk[1]*Mik[2,1]+ Pk[2]*Mik[2,2]+ Pk[3]*Mik[2,3];
  result[3]:= Pk[0]*Mik[3,0]+ Pk[1]*Mik[3,1]+ Pk[2]*Mik[3,2]+ Pk[3]*Mik[3,3];
end;

function mult(const Pk: T_Vect4; const Mki: T_M4): T_Vect4;
var i : shortint;
begin
  for i := 0 to 3 do
  begin
   result[i]:=
     Pk[0]*Mki[0,i]+
     Pk[1]*Mki[1,i]+
     Pk[2]*Mki[2,i]+
     Pk[3]*Mki[3,i];
  end;
end;

function mult(const Mik : T_Tens; const c: real): T_Tens;
begin
  Result.x.x := Mik.x.x*c;
  Result.y.x := Mik.y.x*c;
  Result.z.x := Mik.z.x*c;

  Result.x.y := Mik.x.y*c;
  Result.y.y := Mik.y.y*c;
  Result.z.y := Mik.z.y*c;

  Result.x.z := Mik.x.z*c;
  Result.y.z := Mik.y.z*c;
  Result.z.z := Mik.z.z*c;
end;

function mult(const Mij: T_Tens; const Vj: T_Vect): T_Vect;
begin
  { Ri = Mij*Vj }
  Result.x := Mij.x.x*Vj.x + Mij.x.y*Vj.y+Mij.x.z*Vj.z;
  Result.y := Mij.y.x*Vj.x + Mij.y.y*Vj.y+Mij.y.z*Vj.z;
  Result.z := Mij.z.x*Vj.x + Mij.z.y*Vj.y+Mij.z.z*Vj.z;

end;

function mult_transponed(const Mik: T_Tens; const Mjk : T_Tens): T_Tens;
begin
  // Result.i.j = A.i.k* B.j.k
  Result.x.x := Mik.x.x*Mjk.x.x + Mik.x.y*Mjk.x.y + Mik.x.z*Mjk.x.z;
  Result.x.y := Mik.x.x*Mjk.y.x + Mik.x.y*Mjk.y.y + Mik.x.z*Mjk.y.z;
  Result.x.z := Mik.x.x*Mjk.z.x + Mik.x.y*Mjk.z.y + Mik.x.z*Mjk.z.z;

  Result.y.x := Mik.y.x*Mjk.x.x + Mik.y.y*Mjk.x.y + Mik.y.z*Mjk.x.z;
  Result.y.y := Mik.y.x*Mjk.y.x + Mik.y.y*Mjk.y.y + Mik.y.z*Mjk.y.z;
  Result.y.z := Mik.y.x*Mjk.z.x + Mik.y.y*Mjk.z.y + Mik.y.z*Mjk.z.z;

  Result.z.x := Mik.z.x*Mjk.x.x + Mik.z.y*Mjk.x.y + Mik.z.z*Mjk.x.z;
  Result.z.y := Mik.z.x*Mjk.y.x + Mik.z.y*Mjk.y.y + Mik.z.z*Mjk.y.z;
  Result.z.z := Mik.z.x*Mjk.z.x + Mik.z.y*Mjk.z.y + Mik.z.z*Mjk.z.z;


end;

function mult(const a : T_vect; const Lmd : real): T_Vect;
begin
  result.x := a.x*Lmd;
  result.y := a.y*Lmd;
  result.z := a.z*Lmd;
end;

function mult(const Pk: T_Vect; const Mki : T_M4): T_Vect;
begin
// строка на столбец
// умножает вектор на расширенную матрицу трансформации,
// у которой в последней строке лежат коэффициенты линейного переноса.
   result.x:= Pk.x*Mki[0,0]+ Pk.y*Mki[1,0]+ Pk.z*Mki[2,0]+ Mki[3,0];
   result.y:= Pk.x*Mki[0,1]+ Pk.y*Mki[1,1]+ Pk.z*Mki[2,1]+ Mki[3,1];
   result.z:= Pk.x*Mki[0,2]+ Pk.y*Mki[1,2]+ Pk.z*Mki[2,2]+ Mki[3,2];
end;

function mult(const Vk : T_Vect;const Mki : T_Tens): T_Vect;
begin
// ROW to COLUMN, stroka na stolbec
  Result.x := Vk.x*Mki.x.x+
               Vk.y*Mki.y.x+
               Vk.z*Mki.z.x ;

  Result.y := Vk.x*Mki.x.y+
               Vk.y*Mki.y.y+
               Vk.z*Mki.z.y ;

  Result.z := Vk.x*Mki.x.z+
               Vk.y*Mki.y.z+
               Vk.z*Mki.z.z ;
end;

function mult(const V0,V1 : T_Vect): T_Vect;
begin
  Result.x := V0.y*V1.z- V0.z*V1.y;
  Result.y := V0.z*V1.x- V0.x*V1.z;
  Result.z := V0.x*V1.y- V0.y*V1.x;
end;

function Norma(const DC : T_Vect4): real;
begin
  Result :=  Sqrt ( Sqr(DC[0])  +Sqr(DC[1])  +Sqr(DC[2])  +Sqr(DC[3]) );
end;

function Norma(const S : T_M4): real;
begin
  Result :=  Sqrt ( Sqr(S[0,0])  +Sqr(S[0,1])  +Sqr(S[0,2])  +Sqr(S[0,3]) +
                    Sqr(S[1,0])  +Sqr(S[1,1])  +Sqr(S[1,2])  +Sqr(S[1,3]) +
                    Sqr(S[2,0])  +Sqr(S[2,1])  +Sqr(S[2,2])  +Sqr(S[2,3]) +
                    Sqr(S[3,0])  +Sqr(S[3,1])  +Sqr(S[3,2])  +Sqr(S[3,3]) )
                   ;
end;

function Norma(const S : T_M4; const Diag : real): real;
begin
  Result :=  Sqrt ( Sqr(S[0,0]-Diag)  +Sqr(S[0,1])  +Sqr(S[0,2])  +Sqr(S[0,3]) +
                    Sqr(S[1,0])  +Sqr(S[1,1]-Diag)  +Sqr(S[1,2])  +Sqr(S[1,3]) +
                    Sqr(S[2,0])  +Sqr(S[2,1])  +Sqr(S[2,2]-Diag)  +Sqr(S[2,3]) +
                    Sqr(S[3,0])  +Sqr(S[3,1])  +Sqr(S[3,2])  +Sqr(S[3,3]-Diag) )
                   ;
end;

function Norma(const S1,S2 : T_M4): real;
begin
  Result :=  Sqrt ( Sqr(S1[0 , 0] -S2[ 0,0])  +Sqr(S1[0 , 1] -S2[ 0,1])  +Sqr(S1[0 , 2] -S2[ 0,2])  +Sqr(S1[0 , 3] -S2[ 0,3]) +
                    Sqr(S1[1 , 0] -S2[ 1,0])  +Sqr(S1[1 , 1] -S2[ 1,1])  +Sqr(S1[1 , 2] -S2[ 1,2])  +Sqr(S1[1 , 3] -S2[ 1,3]) +
                    Sqr(S1[2 , 0] -S2[ 2,0])  +Sqr(S1[2 , 1] -S2[ 2,1])  +Sqr(S1[2 , 2] -S2[ 2,2])  +Sqr(S1[2 , 3] -S2[ 2,3]) +
                    Sqr(S1[3 , 0] -S2[ 3,0])  +Sqr(S1[3 , 1] -S2[ 3,1])  +Sqr(S1[3 , 2] -S2[ 3,2])  +Sqr(S1[3 , 3] -S2[ 3,3]) )
                   ;
end;

function Norma(const S : T_Tens): real;
begin
  Result :=  Sqrt ( Sqr(S.x.x)  +Sqr(S.x.y)  +Sqr(S.x.z)  +
                    Sqr(S.y.x)  +Sqr(S.y.y)  +Sqr(S.y.z)  +
                    Sqr(S.z.x)  +Sqr(S.z.y)  +Sqr(S.z.z)  )
                   ;
end;

function Norma(const V : T_Vect): real;
begin
  Result :=  Sqrt ( Sqr(V.x)  +Sqr(V.y)  +Sqr(V.z) ) ;
end;

function Norma(const V1,V2 : T_Vect): real;
begin
  Result :=  Sqrt ( Sqr(V1.x-V2.x)  +
                    Sqr(V1.y-V2.y)  +
                    Sqr(V1.z-V2.z) ) ;
end;

function razn(const a,b : real): real;
begin
  Result := a - b;
end;

function razn(const a,b : T_Tens): T_Tens;
begin
  Result.x.x := a.x.x-b.x.x;
  Result.x.y := a.x.y-b.x.y;
  Result.x.z := a.x.z-b.x.z;

  Result.y.x := a.y.x-b.y.x;
  Result.y.y := a.y.y-b.y.y;
  Result.y.z := a.y.z-b.y.z;

  Result.z.x := a.z.x-b.z.x;
  Result.z.y := a.z.y-b.z.y;
  Result.z.z := a.z.z-b.z.z;

  {$IFDEF ALIGNED_VECTORS_ON}
    Result[0,3] := 0.0 ;
    Result[1,3] := 0.0 ;
    Result[2,3] := 0.0 ;


   {$ENDIF}

end;

function razn(const a,b : T_Vect): T_Vect;
begin
  Result.x := a.x-b.x;
  Result.y := a.y-b.y;
  Result.z := a.z-b.z;
  {$IFDEF ALIGNED_VECTORS_ON} Result[3] := 0.0 ; {$ENDIF}
end;

procedure SetUnite(out a : T_M4);
begin
  a[0,  0] := 1.0;  a[0,  1] := 0.0;  a[0,  2] := 0.0;  a[0,  3] := 0.0;
  a[1,  0] := 0.0;  a[1,  1] := 1.0;  a[1,  2] := 0.0;  a[1,  3] := 0.0;
  a[2,  0] := 0.0;  a[2,  1] := 0.0;  a[2,  2] := 1.0;  a[2,  3] := 0.0;
  a[3,  0] := 0.0;  a[3,  1] := 0.0;  a[3,  2] := 0.0;  a[3,  3] := 1.0;


end;

procedure SetUnite(out a: T_Tens);
begin
  a.x.x := 1.0;
  a.x.y := 0.0;
  a.x.z := 0.0;

  a.y.x := 0.0;
  a.y.y := 1.0;
  a.y.z := 0.0;

  a.z.x := 0.0;
  a.z.y := 0.0;
  a.z.z := 1.0;

end;

procedure set_bit(var b: byte; const bit: byte; const st: boolean);
begin
  if st
  then
    b := b or bit
  else
    b := b and (not bit);


end;

function summ(const a: array of T_Vect): T_Vect;
var i , n : integer;
begin
  n := high(a);
  result := a[0];
  for i := 1 to n do
  begin
    Result.x := Result.x +a[i].x;
    Result.y := Result.y +a[i].y;
    Result.z := Result.z +a[i].z;
  end;
end;

function summ(const a: T_Tens; const b : real): T_Tens;
begin
  Result.x.x := a.x.x+b;
  Result.y.x := a.y.x;
  Result.z.x := a.z.x;

  Result.x.y := a.x.y;
  Result.y.y := a.y.y+b;
  Result.z.y := a.z.y;

  Result.x.z := a.x.z;
  Result.y.z := a.y.z;
  Result.z.z := a.z.z+b;
end;

function summ(const a,b : T_Tens): T_Tens;
begin
  Result.x.x := a.x.x+b.x.x;
  Result.y.x := a.y.x+b.y.x;
  Result.z.x := a.z.x+b.z.x;

  Result.x.y := a.x.y+b.x.y;
  Result.y.y := a.y.y+b.y.y;
  Result.z.y := a.z.y+b.z.y;

  Result.x.z := a.x.z+b.x.z;
  Result.y.z := a.y.z+b.y.z;
  Result.z.z := a.z.z+b.z.z;
end;

function summ(const a: T_Vect; const b : real): T_Vect;
begin
  Result.x := a.x+b;
  Result.y := a.y+b;
  Result.z := a.z+b;
end;

function summ(const a,b : T_Vect): T_Vect;
begin
  Result.x := a.x+b.x;
  Result.y := a.y+b.y;
  Result.z := a.z+b.z;
end;

function summ_t(const a: array of T_Tens): T_Tens;
var i , n : integer;
begin
  n := high(a);
  result := a[0];
  for i := 1 to n do
  begin
    Result.x.x := Result.x.x +a[i].x.x;
    Result.y.x := Result.y.x +a[i].y.x;
    Result.z.x := Result.z.x +a[i].z.x;

    Result.x.y := Result.x.y +a[i].x.y;
    Result.y.y := Result.y.y +a[i].y.y;
    Result.z.y := Result.z.y +a[i].z.y;

    Result.x.z := Result.x.z +a[i].x.z;
    Result.y.z := Result.y.z +a[i].y.z;
    Result.z.z := Result.z.z +a[i].z.z;
  end;
end;

procedure swap(var a,b: real);
var t: real;
begin
  t := b;
  b := a;
  a := t;
{  PInt64(@a)^ := PInt64(@a)^ xor PInt64(@b)^;
  PInt64(@b)^ := PInt64(@a)^ xor PInt64(@b)^;
  PInt64(@a)^ := PInt64(@a)^ xor PInt64(@b)^;
{  a := a -b; //  a-1 b1
  b := b +a; //  a-1 b0
  a := b -a;  //  a1 b0}
end;

function svertka(const M0, M1 : T_Tens): real;
begin
  result := M0.x.x*M1.x.x+M0.x.y*M1.x.y+M0.x.z*M1.x.z+
            M0.y.x*M1.y.x+M0.y.y*M1.y.y+M0.y.z*M1.y.z+
            M0.z.x*M1.z.x+M0.z.y*M1.z.y+M0.z.z*M1.z.z
            ;
end;

function tensMax(const a,b : T_Tens): T_Tens;
begin
  Result.x.x :=  Max(a.x.x, b.x.x);
  Result.x.y :=  Max(a.x.y, b.x.y);
  Result.x.z :=  Max(a.x.z, b.x.z);

  Result.y.x :=  Max(a.y.x, b.y.x);
  Result.y.y :=  Max(a.y.y, b.y.y);
  Result.y.z :=  Max(a.y.z, b.y.z);

  Result.z.x :=  Max(a.z.x, b.z.x);
  Result.z.y :=  Max(a.z.y, b.z.y);
  Result.z.z :=  Max(a.z.z, b.z.z);
end;

function Tensor(const DiagonalVal : real): T_Tens;
begin
  result.x.x := DiagonalVal;
  result.x.y := 0.0;
  result.x.z := 0.0;


  result.y.x := 0.0;
  result.y.y := DiagonalVal;
  result.y.z := 0.0;

  result.z.x := 0.0;
  result.z.y := 0.0;
  result.z.z := DiagonalVal;

  result.x.t := 0.0;
  result.y.t := 0.0;
  result.z.t := 0.0;

end;

function Tensor(const t00,t01,t02, t10,t11,t12, t20,t21,t22 : real): T_Tens;
begin
  result .x.x := t00;
  result .x.y := t01;
  result .x.z := t02;

  result .y.x := t10;
  result .y.y := t11;
  result .y.z := t12;

  result .z.x := t20;
  result .z.y := t21;
  result .z.z := t22;

  result.x.t := 0.0;
  result.y.t := 0.0;
  result.z.t := 0.0;
end;

function Tensor(const F0,F1,F2 : T_Vect): T_Tens;
begin
  result.x := F0;
  result.y := F1;
  result.z := F2;
end;

procedure Transpone(var S: T_M4);
var t: real;
begin
  t := S[0,1]; S[0,1] := S[1,0];  S[1,0] := t;
  t := S[0,2]; S[0,2] := S[2,0];  S[2,0] := t;
  t := S[0,3]; S[0,3] := S[3,0];  S[3,0] := t;

  t := S[1,2]; S[1,2] := S[2,1];  S[2,1] := t;
  t := S[1,3]; S[1,3] := S[3,1];  S[3,1] := t;

  t := S[2,3]; S[2,3] := S[3,2];  S[3,2] := t;


end;

procedure Transpone(var S: T_Tens);
var t: real;
begin
  t := S.x.y; S.x.y := S.y.x;  S.y.x := t;
  t := S.x.z; S.x.z := S.z.x;  S.z.x := t;
  t := S.z.y; S.z.y := S.y.z;  S.y.z := t;
end;

function Vect(const X,Y,Z : real): T_Vect;
begin
  result.x := X;
  result.y := Y;
  result.z := Z;
  result.t := 0.0;
 {$IFDEF ALIGNED_VECTORS_ON}
  result[3] := 0.0;
  {$ENDIF}
end;

function Vect(const DC : T_Vect4): T_Vect;
begin
  result.x := DC[0];
  result.y := DC[1];
  result.z := DC[2];
 {$IFDEF ALIGNED_VECTORS_ON}
  result[3] := 0.0;
  {$ENDIF}
end;

function VectCube(const a : T_Vect): real;
var t1 : real;
begin
  t1 :=sqr(a.x)+sqr(a.y)+sqr(a.z);
  result := sqrt(t1)*t1;
end;

function VectCube(const a,b : T_Vect): real;
var t1 : real;
begin
  t1 :=sqr(a.x-b.x)+sqr(a.y-b.y)+sqr(a.z-b.z);
  result := sqrt(t1)*t1;


end;

function Vector4(const v0,v1,v2,v3: real): T_Vect4;
begin
  T_Vect4(result)[0] := v0;
  T_Vect4(result)[1] := v1;
  T_Vect4(result)[2] := v2;
  T_Vect4(result)[3] := v3;
end;

function Volume(const S: T_Tens): real;
// returns Volume of triangle piramid = 1/6 parallelepiped volume!
begin

  Result :=( S.x.x*(S.y.y*S.z.z - S.z.y*S.y.z )
            -S.y.x*(S.x.y*S.z.z - S.z.y*S.x.z )
            +S.z.x*(S.x.y*S.y.z - S.y.y*S.x.z ) ) / 6;
end;

function Volume(const S0,S1,S2: T_Vect): real;
begin
  Result := ( S0.x*(S1.y*S2.z - S2.y*S1.z )
             -S1.x*(S0.y*S2.z - S2.y*S0.z )
             +S2.x*(S0.y*S1.z - S1.y*S0.z ) )/6 ;
end;

function Volume(const center, S0,S1,S2: T_Vect): real;
begin
  Result := Volume(S0,S1,S2) -
            Volume(center,S1,S2) -
            Volume(S0,center,S2) -
            Volume(S0,S1,center);
end;

function VSqr(const M4: T_M4): real;
begin
  Result := sqr( M4[0,0] )+sqr( M4[0,1] )+sqr( M4[0,2] )+sqr( M4[0,3] )+
            sqr( M4[1,0] )+sqr( M4[1,1] )+sqr( M4[1,2] )+sqr( M4[1,3] )+
            sqr( M4[2,0] )+sqr( M4[2,1] )+sqr( M4[2,2] )+sqr( M4[2,3] )+
            sqr( M4[3,0] )+sqr( M4[3,1] )+sqr( M4[3,2] )+sqr( M4[3,3] );
end;

function VSqr(const a,b : T_Tens): real;
begin
  Result := VSqr( Razn(a, b) );
end;

function VSqr(const a: T_Vect): real;
begin
  Result := sqr( a.x ) + sqr( a.y ) + sqr( a.z );
end;

function VSqr(const a,b: T_Vect): real;
begin
  Result := sqr( a.x-b.x ) + sqr( a.y-b.y ) + sqr( a.z-b.z );
end;

function DotProd(const V0,V1 : T_Vect4): real;
begin result := V0[0] * V1[0] + V0[1] * V1[1]+ V0[2] * V1[2]+ V0[3] * V1[3];
end;



procedure Zero(out a : integer);
begin
  a := 0;
end;

procedure Zero(out a : real);
begin
  a := 0.0;
end;

procedure Zero(out a : T_M4);
begin
  a[ 0,  0] := 0.0;  a[ 0,  1] := 0.0;  a[ 0,  2] := 0.0;  a[ 0, 3] := 0.0;
  a[ 1,  0] := 0.0;  a[ 1,  1] := 0.0;  a[ 1,  2] := 0.0;  a[ 1, 3] := 0.0;
  a[ 2,  0] := 0.0;  a[ 2,  1] := 0.0;  a[ 2,  2] := 0.0;  a[ 2, 3] := 0.0;
  a[ 3,  0] := 0.0;  a[ 3,  1] := 0.0;  a[ 3,  2] := 0.0;  a[ 3, 3] := 0.0;
end;

procedure Zero(out a : T_Vect4);
begin
  a[ 0] := 0.0;
  a[ 1] := 0.0;
  a[ 2] := 0.0;
  a[ 3] := 0.0;
end;

procedure Zero(out a : T_RealArr);
var i: integer;
begin
  for i := 0 to High(a) do
    a[i] := 0.0;


end;

procedure Zero(out a : T_Tens);
begin
  a.x.x := 0.0;   a.x.y := 0.0;   a.x.z := 0.0;
  a.y.x := 0.0;   a.y.y := 0.0;   a.y.z := 0.0;
  a.z.x := 0.0;   a.z.y := 0.0;   a.z.z := 0.0;
end;

procedure Zero(out a : T_TensorArr);
var i: integer;
begin
  for i := 0 to High(a) do
    ZERO(a[i]);
end;

procedure Zero(out a : T_Vect);
begin
  a.x := 0.0;  a.y := 0.0;  a.z := 0.0;  a.t := 0.0;
end;

procedure Zero(out a : T_VectArr);
var i: integer;
begin
  for i := 0 to High(a) do
  begin
    a[i].x := 0.0;
    a[i].y := 0.0;
    a[i].z := 0.0;
    a[i].t := 0.0;
  end;


end;

procedure Add(var B_ab: T_Tens; const _DR_alfa, _DR_beta: T_Vect; const alfa:
    real);
begin
    // B[alfa, beta ] = _DR[alfa]*DC_beta,  dR[-1] = 1.0;

    B_ab.x.x := B_ab.x.x + _DR_alfa.x * _DR_beta.x*alfa;
    B_ab.x.y := B_ab.x.y + _DR_alfa.x * _DR_beta.y*alfa;
    B_ab.x.z := B_ab.x.z + _DR_alfa.x * _DR_beta.z*alfa;

    B_ab.y.x := B_ab.y.x + _DR_alfa.y * _DR_beta.x*alfa;
    B_ab.y.y := B_ab.y.y + _DR_alfa.y * _DR_beta.y*alfa;
    B_ab.y.z := B_ab.y.z + _DR_alfa.y * _DR_beta.z*alfa;

    B_ab.z.x := B_ab.z.x + _DR_alfa.z * _DR_beta.x*alfa;
    B_ab.z.y := B_ab.z.y + _DR_alfa.z * _DR_beta.y*alfa;
    B_ab.z.z := B_ab.z.z + _DR_alfa.z * _DR_beta.z*alfa;
end;

function Norma(const V1,V2 : T_RealArr; const cnt: integer = -1): real;
var i: integer;
    count : integer;
begin
  if cnt >=0 then count := cnt
  else count := Length(V1);

  result := 0;

  for i := 0 to count -1 do
    result := result + sqr(V1[i]-V2[i]);

  result := sqrt(result);

end;

procedure MaxVect(var a: T_Vect; const b: T_Vect);
begin
  if a.x < b.x then a.x := b.x;
  if a.y < b.y then a.y := b.y;
  if a.z < b.z then a.z := b.z;
end;

procedure MinVect(var a: T_Vect; const b: T_Vect);
begin
  if a.x > b.x then a.x := b.x;
  if a.y > b.y then a.y := b.y;
  if a.z > b.z then a.z := b.z;
end;

procedure Invert(var S: T_Tens);
var  D : real;
     Src: T_Tens;
begin

  Src := S;
  D :=(   Src.x.x*( Src.y.y*Src.z.z - Src.y.z*Src.z.y )
        + Src.x.y*( Src.z.x*Src.y.z - Src.y.x*Src.z.z )
        + Src.x.z*( Src.y.x*Src.z.y - Src.y.y*Src.z.x ) );

  if Abs(d) < 1.0e-10 then
  begin
    Exit;
  end;


  D:= 1.0/ D;

  S.x.x :=    (  Src.y.y*Src.z.z - Src.y.z*Src.z.y  )* D;
  S.x.y :=   -(  Src.x.y*Src.z.z - Src.x.z*Src.z.y  )* D;
  S.x.z :=   +(  Src.x.y*Src.y.z - Src.x.z*Src.y.y  )* D;

  S.y.x :=   -(  Src.y.x*Src.z.z - Src.z.x*Src.y.z  )* D;
  S.y.y :=   +(  Src.x.x*Src.z.z - Src.x.z*Src.z.x  )* D;
  S.y.z :=   -(  Src.x.x*Src.y.z - Src.y.x*Src.x.z  )* D;

  S.z.x :=    (  Src.y.x*Src.z.y - Src.y.y*Src.z.x  )* D;
  S.z.y :=   -(  Src.x.x*Src.z.y - Src.x.y*Src.z.x  )* D;
  S.z.z :=   +(  Src.x.x*Src.y.y - Src.x.y*Src.y.x  )* D;


end;

function Norma(const V : T_RealArr; const cnt: integer = -1): real;
var i: integer;
    count : integer;
    Res: real;
begin
  if cnt >=0 then count := cnt
  else count := Length(V);

  res := 0;

  for i := 0 to count -1 do
    res := res + sqr(V[i]);

  result := sqrt(res);

end;

procedure Add(var a: T_Vect; const C: real; const b: T_Vect; const d: real);
begin
  a.x := a.x+b.x*C*d;
  a.y := a.y+b.y*C*d;
  a.z := a.z+b.z*C*d;
end;

procedure Add(var a: T_Vect; const C: real; const b: T_Vect);
begin
  a.x := a.x+b.x*C;
  a.y := a.y+b.y*C;
  a.z := a.z+b.z*C;
end;

function mult(const Mij: T_Tens; const Vj, Vi: T_Vect): real;
begin
  { Ri = Mij*Vj }
  Result :=
      (Mij.x.x*Vj.x + Mij.x.y*Vj.y+Mij.x.z*Vj.z)*Vi.x+
      (Mij.y.x*Vj.x + Mij.y.y*Vj.y+Mij.y.z*Vj.z)*Vi.y+
      (Mij.z.x*Vj.x + Mij.z.y*Vj.y+Mij.z.z*Vj.z)*Vi.z;

end;

procedure Add(var a: real; const b, c : T_Vect; const V: real);
begin
  a := a + V * DotProd( b , c );
end;

procedure Add(var B_ab: T_Tens; const _DR_alfa, _DR_beta: T_Vect);
begin
    // B[alfa, beta ] = _DR[alfa]*DC_beta,  dR[-1] = 1.0;

    B_ab.x.x := B_ab.x.x + _DR_alfa.x * _DR_beta.x;
    B_ab.x.y := B_ab.x.y + _DR_alfa.x * _DR_beta.y;
    B_ab.x.z := B_ab.x.z + _DR_alfa.x * _DR_beta.z;


    B_ab.y.x := B_ab.y.x + _DR_alfa.y * _DR_beta.x;
    B_ab.y.y := B_ab.y.y + _DR_alfa.y * _DR_beta.y;
    B_ab.y.z := B_ab.y.z + _DR_alfa.y * _DR_beta.z;


    B_ab.z.x := B_ab.z.x + _DR_alfa.z * _DR_beta.x;
    B_ab.z.y := B_ab.z.y + _DR_alfa.z * _DR_beta.y;
    B_ab.z.z := B_ab.z.z + _DR_alfa.z * _DR_beta.z;

end;

function DotProd(const V0,V1 : T_Vect): real;
begin
  result := V0.x * V1.x + V0.y * V1.y+ V0.z * V1.z;
end;

function sse_array_svertka(const V0: pdouble; const N: integer): real;
    register;
asm

 {$IFDEF CPUX64}
    .NOFRAME

    // RCX- V0
    // RDX- N


    @rep:

     pxor    xmm0, xmm0 // XMM0 = zero
     pxor    xmm1, xmm1
     pxor    xmm2, xmm2
     pxor    xmm3, xmm3
     pxor    xmm4, xmm4


     sub     rdx,  1

     movapd   xmm2, OWORD [rcx ]         // Закидываем первый вектор в регистры
     movapd   xmm3, OWORD [rcx+16 ]

    @next:
     add      rcx,  32
     movapd   xmm0, xmm2           // пропихиваем вектор вниз
     movapd   xmm1, xmm3
     movapd   xmm2, OWORD [rcx ]     // докидываем новый вектор
     movapd   xmm3, OWORD [rcx+16 ]

     mulpd XMM0, XMM2        // перемножаем первую пару компонент
     mulpd XMM1, XMM3        // перемножаем вторую пару компонент

     addpd  XMM4, XMM0       // добавляем в результат первую пару компонент
     addpd  XMM4, XMM1       // -- вторую пару компонент

     dec     rdx             // уменьшаем коунтер
     jnz     @next           // прыг на новую итерацию.


     haddpd  XMM4, XMM4      // горизонтально суммируем полукомпоненты
     movapd  xmm0, xmm4      // result of float calc returned via XMM0!
 {$ELSE}

    @rep:
    {Вычисление суммы}
//     mov     rcx,  [N]
     pxor    xmm0, xmm0 // XMM0 = zero
     pxor    xmm1, xmm1
     pxor    xmm2, xmm2
     pxor    xmm3, xmm3
     pxor    xmm4, xmm4


     sub     edx,  1

     movupd   xmm2, OWORD [eax ]
     movupd   xmm3, OWORD [eax+16 ]

    @next:
     add      eax,  32
     movapd   xmm0, xmm2
     movapd   xmm1, xmm3
     movupd   xmm2, OWORD [eax ]
     movupd   xmm3, OWORD [eax+16 ]

     mulpd XMM0, XMM2
     mulpd XMM1, XMM3

     addpd  XMM4, XMM0
     addpd  XMM4, XMM1

     dec     edx
     jnz     @next


     haddpd  XMM4, XMM4
     movapd  xmm0, xmm4 // result of float calc returned via XMM0!
     movsd   qWORD [Result], xmm0
 {$ENDIF}

end;

procedure Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j: T_Vect);
begin
  a_i.x := a_i.x+b_ij.x.x*C_j.x+b_ij.x.y*C_j.y+b_ij.x.z*C_j.z ;
  a_i.y := a_i.y+b_ij.y.x*C_j.x+b_ij.y.y*C_j.y+b_ij.y.z*C_j.z ;
  a_i.z := a_i.z+b_ij.z.x*C_j.x+b_ij.z.y*C_j.y+b_ij.z.z*C_j.z ;
end;

procedure Invert_gauss(var S: T_Tens);
var  D : real;
begin
// STRAIGTH MOVEMENT!
  // STEP  #1

  d := 1.0 / S.x.x;
  S.x.x :=       D;
  S.x.y := S.x.y*D;
  S.x.z := S.x.z*D;



    // STEP  #2
  D := -S.y.x;
  s.y.x :=         S.x.x*d;
  s.y.y := S.y.y + S.x.y*d;
  s.y.z := S.y.z + S.x.z*d;


    // STEP  #3
  D := -S.z.x;
  S.z.x :=         S.x.x*d;
  S.z.y := S.z.y + S.x.y*d;
  S.z.z := S.z.z + S.x.z*d;


    // STEP  #4
  D := 1.0/S.y.y;
  S.y.x := S.y.x*D;
  S.y.y :=       D;
  S.y.z := S.y.z*D;


    // STEP  #5
  D := -S.z.y;
  S.z.x := S.z.x + S.y.x*D;
  S.z.y :=         S.y.y*D;
  S.z.z := S.z.z + S.y.z*D;


    // STEP  #6
  D := 1.0/S.z.z;
  S.z.x := S.z.x*D;
  S.z.y := S.z.y*D;
  S.z.z :=       D;    // 1/A22

// BACKWARD MOVEMENT!

    // STEP  #7
  D := -S.x.y;
  S.x.x := S.x.x + S.y.x*D;
  S.x.y :=         S.y.y*D;
  S.x.z := S.x.z + S.y.z*D;



    // STEP  #8
  D := -S.y.z;
  S.y.x := S.y.x + S.z.x*D;
  S.y.y := S.y.y + S.z.y*D;
  S.y.z :=         S.z.z*D;


    // STEP  #9
  D := -S.x.z;
  S.x.x := S.x.x + S.z.x*D;
  S.x.y := S.x.y + S.z.y*D;
  S.x.z :=         S.z.z*D;

end;

procedure Invert_delphi_assembler_cleared(var S: T_Tens);
const Sign: UInt64= $7FFFFFFFFFFFFFFF;
      eps: double = 1.0e-10;
asm
push rbp
sub rsp,$0000000000000090
mov rbp,rsp
mov [rbp+$000000a0],rcx
//sph3d.base.types.pas.5188: Res.x.x :=    ( + S.y.y*S.z.z - S.y.z*S.z.y  );
mov rax,[rbp+$000000a0]
mov rcx,[rbp+$000000a0]
movsd xmm0,qword ptr [rax+$28]
mulsd xmm0,qword ptr [rcx+$50]
mov rax,[rbp+$000000a0]
mov rcx,[rbp+$000000a0]
movsd xmm1,qword ptr [rax+$30]
mulsd xmm1,qword ptr [rcx+$48]
subsd xmm0,xmm1
movsd qword ptr [rbp+$28],xmm0
//sph3d.base.types.pas.5189: Res.x.y :=    ( + S.x.z*S.z.y - S.x.y*S.z.z  );
movsd xmm0,qword ptr [rax+$10]
mulsd xmm0,qword ptr [rcx+$48]
movsd xmm1,qword ptr [rax+$08]
mulsd xmm1,qword ptr [rcx+$50]
subsd xmm0,xmm1
movsd qword ptr [rbp+$30],xmm0
//sph3d.base.types.pas.5191: Res.x.z :=    ( + S.x.y*S.y.z - S.x.z*S.y.y  );
movsd xmm0,qword ptr [rax+$08]
mulsd xmm0,qword ptr [rcx+$30]
movsd xmm1,qword ptr [rax+$10]
mulsd xmm1,qword ptr [rcx+$28]
subsd xmm0,xmm1
movsd qword ptr [rbp+$38],xmm0
//sph3d.base.types.pas.5193: Res.y.x :=    ( + S.z.x*S.y.z - S.y.x*S.z.z  );
movsd xmm0,qword ptr [rax+$40]
mulsd xmm0,qword ptr [rcx+$30]
movsd xmm1,qword ptr [rax+$20]
mulsd xmm1,qword ptr [rcx+$50]
subsd xmm0,xmm1
movsd qword ptr [rbp+$48],xmm0
//sph3d.base.types.pas.5194: Res.y.y :=    ( + S.x.x*S.z.z - S.x.z*S.z.x  );
movsd xmm0,qword ptr [rax]
mulsd xmm0,qword ptr [rcx+$50]
movsd xmm1,qword ptr [rax+$10]
mulsd xmm1,qword ptr [rcx+$40]
subsd xmm0,xmm1
movsd qword ptr [rbp+$50],xmm0
//sph3d.base.types.pas.5195: Res.y.z :=    ( + S.y.x*S.x.z - S.x.x*S.y.z  );
movsd xmm0,qword ptr [rax+$20]
mulsd xmm0,qword ptr [rcx+$10]
movsd xmm1,qword ptr [rax]
mulsd xmm1,qword ptr [rcx+$30]
subsd xmm0,xmm1
movsd qword ptr [rbp+$58],xmm0
//sph3d.base.types.pas.5197: Res.z.x :=    ( + S.y.x*S.z.y - S.y.y*S.z.x  );
movsd xmm0,qword ptr [rax+$20]
mulsd xmm0,qword ptr [rcx+$48]
movsd xmm1,qword ptr [rax+$28]
mulsd xmm1,qword ptr [rcx+$40]
subsd xmm0,xmm1
movsd qword ptr [rbp+$68],xmm0
//sph3d.base.types.pas.5198: Res.z.y :=    ( + S.x.y*S.z.x - S.x.x*S.z.y  );
movsd xmm0,qword ptr [rax+$08]
mulsd xmm0,qword ptr [rcx+$40]
movsd xmm1,qword ptr [rax]
mulsd xmm1,qword ptr [rcx+$48]
subsd xmm0,xmm1
movsd qword ptr [rbp+$70],xmm0
//sph3d.base.types.pas.5199: Res.z.z :=    ( + S.x.x*S.y.y - S.x.y*S.y.x  );
movsd xmm0,qword ptr [rax]
mulsd xmm0,qword ptr [rcx+$28]
movsd xmm1,qword ptr [rax+$08]
mulsd xmm1,qword ptr [rcx+$20]
subsd xmm0,xmm1
movsd qword ptr [rbp+$78],xmm0
//sph3d.base.types.pas.5201: D := Res.x.x*S.x.x +

movsd xmm0,qword ptr [rbp+$28]
mulsd xmm0,qword ptr [rax]

movsd xmm1,qword ptr [rbp+$30]
mulsd xmm1,qword ptr [rax+$20]
addsd xmm0,xmm1

movsd xmm1,qword ptr [rbp+$38]
mulsd xmm1,qword ptr [rax+$40]
addsd xmm0,xmm1
movsd qword ptr [rbp+$00000088],xmm0
//sph3d.base.types.pas.5205: if Abs(d) < 1.0e-10 then
movsd xmm0,qword ptr [rbp+$00000088]
//call @system.abs}
{    movsd XMM0, XMM4     // if ABS (D) > eps then PROCEED else EXIT
    movsd XMM1, sign
    pand XMM0, XMM1
    movsd XMM1, eps
    comisd XMM1, XMM0}

//movsd xmm1,qword ptr [rel $00000112]
//comisd xmm1,xmm0
    movsd XMM1, Sign// $7FFFFFFFFFFFFFFF
    pand XMM0, XMM1
    movsd XMM1, eps //qword ptr [rel $00000112] // eps
    comisd XMM1, XMM0
jnbe Invert_delphi_assembler_cleared + $368
//sph3d.base.types.pas.5208: D := 1.0 / D;
movsd xmm0,qword ptr [rel $00000108]
divsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rbp+$00000088],xmm0
//sph3d.base.types.pas.5212: S.x.x := Res.x.x*D;

movsd xmm0,qword ptr [rbp+$28]
mulsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rax],xmm0
//sph3d.base.types.pas.5213: S.x.y := Res.x.y*D;

movsd xmm0,qword ptr [rbp+$30]
mulsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rax+$08],xmm0
//sph3d.base.types.pas.5214: S.x.z := Res.x.z*D;

movsd xmm0,qword ptr [rbp+$38]
mulsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rax+$10],xmm0
//sph3d.base.types.pas.5216: S.y.x := Res.y.x*D;
movsd xmm0,qword ptr [rbp+$48]
mulsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rax+$20],xmm0
//sph3d.base.types.pas.5217: S.y.y := Res.y.y*D;
movsd xmm0,qword ptr [rbp+$50]
mulsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rax+$28],xmm0
//sph3d.base.types.pas.5218: S.y.z := Res.y.z*D;
movsd xmm0,qword ptr [rbp+$58]
mulsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rax+$30],xmm0
//sph3d.base.types.pas.5220: S.z.x := Res.z.x*D;
movsd xmm0,qword ptr [rbp+$68]
mulsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rax+$40],xmm0
//sph3d.base.types.pas.5221: S.z.y := Res.z.y*D;
movsd xmm0,qword ptr [rbp+$70]
mulsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rax+$48],xmm0
//sph3d.base.types.pas.5222: S.z.z := Res.z.z*D;
movsd xmm0,qword ptr [rbp+$78]
mulsd xmm0,qword ptr [rbp+$00000088]
movsd qword ptr [rax+$50],xmm0
//sph3d.base.types.pas.5224: end;
lea rsp,[rbp+$00000090]
pop rbp

end;


procedure Invert(const a : T_M4 ; out E : T_M4);
const n = 4;
var   M : array [0..3,0..7] of real;
   i,j,k : integer;
   S1 : real;
   ErrorCode : integer;

begin
  ErrorCode := 0; // никаких ошибок не было.

  for i := 0 to n-1 do
  for j := 0 to n-1 do
  begin
    M[i,j] := a[i,j];
    M[i,j+n] := 0.0;
    if i=j then
      M[i,j+n] := 1.0;
  end;

{Прямой ход. }
  for i := 0 to 3 do
  begin
    if Abs(M[i,i])>1.0e-14
    then begin
      S1 := 1.0/M[i,i];
      for j := 3 downto i+1 do
        M[j,i] := M[j,i]*S1;
    end
    else begin
      ErrorCode := 1; // на главной оси вылез ноль и ничего тут не попишешь.
      Break; //выходим из цикла главного хода.
    end;

    for j := i+1 to 3 do
    begin
      S1 := M[j,i];
      M[j,i] := 0.0;
      for k := i+1 to 7 do
        M[j,k] := M[j,k]-M[i,k]*S1 ;
    end;
  end;

  {Обратный ход}
  if ErrorCode = 0 // если нет ошибoк- то делаем обратный ход
  then begin
    for i := 3 downto 0 do
    begin
      if Abs(M[i,i])>1.0e-14
      then begin
        S1 := 1.0/M[i,i];
        for k := i to 7 do
          M[i,k] := M[i,k]*S1; {Сделали единичку на диагонали }
      end
      else begin
        ErrorCode := 3;
        Break;
      end;

      for j := i-1 downto 0 do
      begin
        S1 := M[j,i];
        for k := i to 7 do
          M[j,k] := M[j,k] - M[i,k]*S1;
      end;
    end;
  end;

  if ErrorCode = 0
  then begin
    for i := 0 to 3 do
    for j := 0 to 3 do
      E[i,j] := M[i,j+n];
  end
  else begin
    SetUnite(T_M4(E));
  end;

end;

procedure Invert(var a : T_M4);
begin
  Invert(a, a);
end;

function Vector3(const v0,v1,v2: real): T_Vect;
begin
  result.x := v0;
  result.y := v1;
  result.z := v2;
end;

function Norma(const S1,S2 : T_Tens): real;
begin
  Result :=  Sqrt ( Sqr(S1.x.x-S2.x.x)  +Sqr(S1.x.y-S2.x.y)  +Sqr(S1.x.z-S2.x.z)  +
                    Sqr(S1.y.x-S2.y.x)  +Sqr(S1.y.y-S2.y.y)  +Sqr(S1.y.z-S2.y.z)  +
                    Sqr(S1.z.x-S2.z.x)  +Sqr(S1.z.y-S2.z.y)  +Sqr(S1.z.z-S2.z.z)  )
                   ;
end;

function VSqr(const a : T_Tens): real;
begin
  Result := sqr( a.x.x ) +  sqr( a.x.y ) +  sqr( a.x.z ) +
            sqr( a.y.x ) +  sqr( a.y.y ) +  sqr( a.y.z ) +
            sqr( a.z.x ) +  sqr( a.z.y ) +  sqr( a.z.z );

end;

procedure Invert_m4(var a: T_M4);
begin
  Invert(a);
end;

procedure Invert2(var S: T_Tens);
var  D : real;
     Res: T_Tens;
begin

  Res.x.x :=    ( + S.y.y*S.z.z - S.y.z*S.z.y  );
  Res.x.y :=    ( + S.x.z*S.z.y - S.x.y*S.z.z  );

  Res.x.z :=    ( + S.x.y*S.y.z - S.x.z*S.y.y  );

  Res.y.x :=    ( + S.z.x*S.y.z - S.y.x*S.z.z  );
  Res.y.y :=    ( + S.x.x*S.z.z - S.x.z*S.z.x  );
  Res.y.z :=    ( + S.y.x*S.x.z - S.x.x*S.y.z  );

  Res.z.x :=    ( + S.y.x*S.z.y - S.y.y*S.z.x  );
  Res.z.y :=    ( + S.x.y*S.z.x - S.x.x*S.z.y  );
  Res.z.z :=    ( + S.x.x*S.y.y - S.x.y*S.y.x  );

  D := Res.x.x*S.x.x +
       Res.x.y*S.y.x +
       Res.x.z*S.z.x;

  if Abs(d) < 1.0e-10 then
    Exit;

  D := 1.0 / D;

//  S := mult(Res,D);

  S.x.x := Res.x.x*D;
  S.x.y := Res.x.y*D;
  S.x.z := Res.x.z*D;

  S.y.x := Res.y.x*D;
  S.y.y := Res.y.y*D;
  S.y.z := Res.y.z*D;

  S.z.x := Res.z.x*D;
  S.z.y := Res.z.y*D;
  S.z.z := Res.z.z*D;

end;

procedure Invert(var S: T_Tens3);
var  D : real;
     Src: T_Tens3;
begin

  Src := S;
  D :=(   Src.xx*( Src.yy*Src.zz - Src.yz*Src.zy )
        + Src.xy*( Src.zx*Src.yz - Src.yx*Src.zz )
        + Src.xz*( Src.yx*Src.zy - Src.yy*Src.zx ) );

  if Abs(d) < 1.0e-10 then
  begin
    Exit;
  end;

  D:= 1.0/ D;

  S.xx :=    (  Src.yy*Src.zz - Src.yz*Src.zy  )* D;
  S.xy :=   -(  Src.xy*Src.zz - Src.xz*Src.zy  )* D;
  S.xz :=   +(  Src.xy*Src.yz - Src.xz*Src.yy  )* D;

  S.yx :=   -(  Src.yx*Src.zz - Src.zx*Src.yz  )* D;
  S.yy :=   +(  Src.xx*Src.zz - Src.xz*Src.zx  )* D;
  S.yz :=   -(  Src.xx*Src.yz - Src.yx*Src.xz  )* D;

  S.zx :=    (  Src.yx*Src.zy - Src.yy*Src.zx  )* D;
  S.zy :=   -(  Src.xx*Src.zy - Src.xy*Src.zx  )* D;
  S.zz :=   +(  Src.xx*Src.yy - Src.xy*Src.yx  )* D;

end;

function mult(const Mik, Mkj : T_Tens): T_Tens;
begin
  Result.x.x := Mik.x.x*Mkj.x.x+ Mik.x.y*Mkj.y.x+  Mik.x.z*Mkj.z.x;
  Result.y.x := Mik.y.x*Mkj.x.x+ Mik.y.y*Mkj.y.x+  Mik.y.z*Mkj.z.x;
  Result.z.x := Mik.z.x*Mkj.x.x+ Mik.z.y*Mkj.y.x+  Mik.z.z*Mkj.z.x;

  Result.x.y := Mik.x.x*Mkj.x.y+ Mik.x.y*Mkj.y.y+  Mik.x.z*Mkj.z.y;
  Result.y.y := Mik.y.x*Mkj.x.y+ Mik.y.y*Mkj.y.y+  Mik.y.z*Mkj.z.y;
  Result.z.y := Mik.z.x*Mkj.x.y+ Mik.z.y*Mkj.y.y+  Mik.z.z*Mkj.z.y;

  Result.x.z := Mik.x.x*Mkj.x.z+ Mik.x.y*Mkj.y.z+  Mik.x.z*Mkj.z.z;
  Result.y.z := Mik.y.x*Mkj.x.z+ Mik.y.y*Mkj.y.z+  Mik.y.z*Mkj.z.z;
  Result.z.z := Mik.z.x*Mkj.x.z+ Mik.z.y*Mkj.y.z+  Mik.z.z*Mkj.z.z;
end;

class operator T_Vect.Subtract(const A, B: T_Vect): T_Vect;
begin
  Result.X := A.x - B.X;
  Result.Y := A.Y - B.Y;
  Result.Z := A.Z - B.Z;
  Result.t := 0;
end;

class operator T_Vect.Add(const A, B: T_Vect): T_Vect;
begin
  Result.X := A.x + B.X;
  Result.Y := A.Y + B.Y;
  Result.Z := A.Z + B.Z;
  Result.t := 0;
end;

class operator T_Vect.Divide(const A: T_Vect; const B: real): T_Vect;
begin
  Result.x := A.x/B;
  Result.y := A.y/B;
  Result.z := A.z/B;
  Result.t := 0;
end;

class operator T_Vect.Multiply(const A, B: T_Vect): real;
begin
  Result := A.x*B.x + A.y*B.y + A.z*B.z;
  // TODO -cMM: T_Vect.Multiply default body inserted
end;

class operator T_Vect.Multiply(const A: T_Vect; const b: real): T_Vect;
begin
  Result.x := A.x*B;
  Result.y := A.y*B;
  Result.z := A.z*B;
  Result.t := 0;
end;

class operator T_Vect.Multiply(const A: real; const b: T_Vect): T_Vect;
begin
  Result.x := B.x*A;
  Result.y := B.y*A;
  Result.z := B.z*A;
  Result.t := 0;
end;

class operator T_Vect.Negative(const A: T_Vect): T_Vect;
begin
  Result.X :=  -A.X;
  Result.Y :=  -A.Y;
  Result.Z :=  -A.Z;
  Result.t := 0;
end;

class operator T_Tens.Multiply(const A, B: T_Tens): T_Tens;
begin
  Result := mult(A, B);

end;

class operator T_Tens.Add(const A: T_Tens; const B: T_Tens): T_Tens;
begin
  Result := Summ(A, B);
end;

class operator T_Tens.Multiply(const A: T_Tens; const B: T_Vect): T_Vect;
begin
  Result := mult(A, B);
end;

class operator T_Tens.Multiply(const A: T_Vect; const B: T_Tens): T_Vect;
begin
  Result := mult(A, B);
end;

class procedure T_SSE.Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j:
    T_Vect; const N: integer);
 {$IFDEF CPUX64}
asm
{  a_i.x := a_i.x+b_ij.x.x*C_j.x+b_ij.x.y*C_j.y+b_ij.x.z*C_j.z ;
  a_i.y := a_i.y+b_ij.y.x*C_j.x+b_ij.y.y*C_j.y+b_ij.y.z*C_j.z ;
  a_i.z := a_i.z+b_ij.z.x*C_j.x+b_ij.z.y*C_j.y+b_ij.z.z*C_j.z ;}
  // C_j: RAX,  R8
  // b_ij:  RDX, r8   (?)
  // a_i: RCX

    .NOFRAME

    .PUSHNV RBX
    .SAVENV XMM5
    .SAVENV XMM6
    .SAVENV XMM7


    XOR RBX, RBX

    mov EBX, N

    @next:
    prefetch  [RDX]
    prefetch  [RDX+2*FloatSize]
    prefetch  [RDX+4*FloatSize]
    prefetch  [RDX+6*FloatSize]
    prefetch  [RDX+8*FloatSize]
    prefetch  [RDX+10*FloatSize]
    prefetch  [RDX+12*FloatSize]
    prefetch  [RCX]
    prefetch  [RCX+2*FloatSize]
    prefetch  [R8]
    prefetch  [R8+2*FloatSize]

    movupd  XMM4, oWORD[a_i]        // XMM4 = a0 | a1
    movupd  XMM5, oWORD[a_i+16]     // XMM5 = a2 | w

    movupd  XMM6, oWORD[C_j]        // XMM6 = c0 | c1
    movupd  XMM7, oWORD[C_j+16]     // XMM7 = c2 | w

    movupd  XMM0, oWORD[ b_ij+0  ]  //  XMM0 = m00 | m01
    movupd  XMM1, oWORD[ b_ij+16 ]  //  XMM1 = m02 | w

    mulpd   XMM0, XMM6              // xmm0 = m00*c0 | m01*c1
    mulsd   XMM1, XMM7              // xmm1 = m02*c2 | 0
    addsd   XMM0, XMM1              // xmm0 = m00*c0 + m01*c1 | m02*c2


    movupd  XMM2, oWORD[ b_ij+32 ]  //  XMM2 = m10 | m11
    movlpd   XMM3, qWORD[ b_ij+48 ]  //  XMM3 = m12 | w

    mulpd   XMM2, XMM6              // xmm2 = m10*c0 | m11*c1
    mulsd   XMM3, XMM7              // xmm3 = m12*c2 | 0
    addpd   XMM2, XMM3              // xmm2 = m10*c0 + m11*c1 | m12*c2

    haddpd  XMM0, XMM2              // xmm0 = m00*c0 + m01*c1 + m02*c2 | m10*c0 + m11*c1 + m12*c2
    addpd   XMM4, XMM0


    movupd  XMM2, oWORD[ b_ij+64 ]  //  XMM2 = m20 | m21
    movlpd  XMM3, qWORD[ b_ij+80 ]  //  XMM3 = m22 | w

    mulpd   XMM2, XMM6              // XMM2 = m20*c0 | m21*c1
    mulsd   XMM3, XMM7              // XMM3 = m22*c2 | 0
    movlhps XMM5, XMM3              // XMM5 = a2, m22*c2

    addpd  XMM5, XMM2              // xmm3 = m22*c2+a2 | m20*c0 + m21*c1
    haddpd  XMM5, XMM5

    movapd oWORD[a_i], XMM4
    movlpd qWORD[a_i+16], XMM5

    add C_j, 4*FloatSize
    add b_ij, 12*FloatSize
    add a_i, 4*FloatSize

    dec     EBX             // уменьшаем коунтер , значение падает в аккумулятор!
    jnz     @next          // прыг на новую итерацию, если аккумулятор <>0




{$ELSE}

var i : integer;
begin

  for i := 0 to N-1 do
  begin
    T_VectArr(@a_i)[i].x := T_VectArr(@a_i)[i].x +
       T_VectArr(@C_j)[i].x * T_TensorArr(@b_ij)[i].x.x +
       T_VectArr(@C_j)[i].y * T_TensorArr(@b_ij)[i].x.y +
       T_VectArr(@C_j)[i].z * T_TensorArr(@b_ij)[i].x.z
       ;

    T_VectArr(@a_i)[i].y := T_VectArr(@a_i)[i].y +
       T_VectArr(@C_j)[i].x * T_TensorArr(@b_ij)[i].y.x +
       T_VectArr(@C_j)[i].y * T_TensorArr(@b_ij)[i].y.y +
       T_VectArr(@C_j)[i].z * T_TensorArr(@b_ij)[i].y.z
       ;

    T_VectArr(@a_i)[i].z := T_VectArr(@a_i)[i].z +
       T_VectArr(@C_j)[i].x * T_TensorArr(@b_ij)[i].z.x +
       T_VectArr(@C_j)[i].y * T_TensorArr(@b_ij)[i].z.y +
       T_VectArr(@C_j)[i].z * T_TensorArr(@b_ij)[i].z.z
       ;

  end;


{$ENDIF CPUX64}

end;

class procedure T_SSE.Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j:
    T_Vect);
{$IFDEF CPUX64}
asm
{  a_i.x := a_i.x+b_ij.x.x*C_j.x+b_ij.x.y*C_j.y+b_ij.x.z*C_j.z ;
  a_i.y := a_i.y+b_ij.y.x*C_j.x+b_ij.y.y*C_j.y+b_ij.y.z*C_j.z ;
  a_i.z := a_i.z+b_ij.z.x*C_j.x+b_ij.z.y*C_j.y+b_ij.z.z*C_j.z ;}
  // C_j: RAX,
  // b_ij:  RDX, r8   (?)
  // a_i: RCX


    .NOFRAME
    .SAVENV XMM5
    .SAVENV XMM6
    .SAVENV XMM7

    movupd  XMM6, oWORD[C_j]         // XMM6 = c0 | c1
    movsd   XMM7, qWORD[C_j+16]      // XMM7 = c2 | w

    movupd  XMM4, oWORD[a_i]         // XMM4 = a0 | a1
    movsd   XMM5, qWORD[a_i+16]      // XMM5 = a2 | w

    movupd  XMM0, oWORD[ b_ij+0  ]  //  XMM0 = m00 | m01
    movsd   XMM1, qWORD[ b_ij+16 ]  //  XMM1 = m02 | w

    movupd  XMM2, oWORD[ b_ij+32 ]  //  XMM2 = m10 | m11
    movsd   XMM3, qWORD[ b_ij+48 ]  //  XMM3 = m12 | w


    mulpd   XMM0, XMM6            // xmm0 = m00*c0 | m01*c1
    mulsd   XMM1, XMM7            // xmm1 = m02*c2 | 0
    haddpd  XMM0, XMM1           // xmm0 = m00*c0 + m01*c1 | m02*c2

    mulpd   XMM2, XMM6            // xmm2 = m10*c0 | m11*c1
    mulsd   XMM3, XMM7            // xmm3 = m22*c2 | 0
    haddpd  XMM2, XMM3           // xmm2 = m10*c0 + m11*c1 | m12*c2


    haddpd  XMM0, XMM2           // xmm0 = m00*c0 + m01*c1 + m02*c2 | m10*c0 + m11*c1 + m12*c2
    addpd   XMM4, XMM0

    movupd  XMM2, oWORD[ b_ij+64 ]  //  XMM2 = m20 | m21
    movsd   XMM3, qWORD[ b_ij+80 ]  //  XMM3 = m22 | w

    mulpd   XMM2, XMM6            // XMM2 = m20*c0 | m21*c1
    mulsd   XMM3, XMM7            // XMM3 = m22*c2 | 0
    movlhps XMM5, XMM3          // XMM5 = a2, m22*c2

    haddpd  XMM5, XMM2           // xmm3 = m22*c2+a2 | m20*c0 + m21*c1
    haddpd  XMM5, XMM5

    movupd  oWORD[a_i], XMM4
    movsd   qWORD[a_i+16], XMM5
{$ELSE}
begin
  a_i.x := a_i.x+b_ij.x.x*C_j.x+b_ij.x.y*C_j.y+b_ij.x.z*C_j.z;
  a_i.y := a_i.y+b_ij.y.x*C_j.x+b_ij.y.y*C_j.y+b_ij.y.z*C_j.z;
  a_i.z := a_i.z+b_ij.z.x*C_j.x+b_ij.z.y*C_j.y+b_ij.z.z*C_j.z;
{$IFEND CPUX64}

end;

class procedure T_SSE.Add(var a_j: T_Vect; const C_i: T_Vect; const b_ij:
    T_Tens);
 {$IFDEF CPUX64}
asm
{  a_i.x := a_i.x+b_ij.x.x*C_j.x+b_ij.x.y*C_j.y+b_ij.x.z*C_j.z ;
  a_i.y := a_i.y+b_ij.y.x*C_j.x+b_ij.y.y*C_j.y+b_ij.y.z*C_j.z ;
  a_i.z := a_i.z+b_ij.z.x*C_j.x+b_ij.z.y*C_j.y+b_ij.z.z*C_j.z ;}
  // C_i: RDX,
  // b_ij:  r8   (?)
  // a_j: RCX

    .NOFRAME


    movupd  XMM4, oWORD[ a_j     ]   // XMM4 = a0 | a1
    movsd   XMM5, qWORD[ a_j+16  ]   // XMM5 = a2 | w

    movapd  XMM0, oWORD[ b_ij  ]     //  XMM0 = m00 | m01
    ADD b_ij, 2*FloatSize            //  R8 = @M+2=  @M.02
    movsd   XMM1, qWORD[ b_ij  ]     //  XMM1 = m02 | w
    ADD b_ij, 2*FloatSize            //  R8 = @M+4   @M.10

    movsd   XMM7, qWORD[ C_i  ]      //  XMM7 = c0 | 0
    ADD     C_i, FloatSize           //  RDX = C_i + 1 elem
    movlhps XMM7, XMM7               //  XMM7 = c0 | c0

    mulpd XMM0, XMM7                 //  xmm0 = m00*c0 | m01*c0
    mulsd XMM1, XMM7                 //  xmm1 = m02*c0 | 0

    addpd  XMM4, XMM0                //  a0+c0*m00 | a1 + c0*m01
    addsd  XMM5, XMM1                //  a2+c0*m02 | 0

    movapd  XMM0, oWORD[ b_ij  ]     //  XMM0 = m10 | m11
    ADD b_ij, 2*FloatSize            //  R8 = @M+6 =@M.12
    movsd   XMM1, qWORD[ b_ij  ]     //  XMM1 = m12 | w
    ADD b_ij, 2*FloatSize            //  R8 = @M+8 =@M.20

    movsd   XMM7, qWORD[ C_i    ]    // XMM7 = c1 | 0
    ADD     C_i, FloatSize           //  RDX = C_i + 2 elem

    movlhps XMM7, XMM7               // XMM7 = c1 | c1

    mulpd XMM0, XMM7                 // xmm0 = m10*c1 | m11*c1
    mulsd XMM1, XMM7                 // xmm1 = m12*c1 | 0

    addpd  XMM4, XMM0                //  a0+c0*m00+m10*c1 | a1 + c0*m01+m11*c1
    addsd  XMM5, XMM1                //  a2+c0*m02+m12*c1 | 0

    movapd  XMM0, oWORD[ b_ij  ]     //  XMM0 = m20 | m21
    ADD b_ij, 2*FloatSize            // R8 = @M+10 = @M.22
    movsd   XMM1, qWORD[ b_ij ]      //  XMM1 = m22 | w
//    Add R8, 16                     // R8 = @M+12 = @M.next ! (12= 3 Vect*4 elems)

    movsd   XMM7, qWORD[ C_i   ]     // XMM7 = c2 | 0
//    Add     RDX, 8                 //  RDX = C_i + 3 elem
    movlhps XMM7, XMM7               // XMM7 = c2 | c2

    mulpd XMM0, XMM7                 // xmm0 = m20*c2 | m21*c2
    mulsd XMM1, XMM7                 // xmm1 = m22*c2 | 0

    addpd  XMM4, XMM0                //  a0+c0*m00+m10*c1+m20*c2 | a1 + c0*m01+m11*c1+m21*c2
    addsd  XMM5, XMM1                //  a2+c0*m02+m12*c1+m22*c2 | 0

    movupd  oWORD[a_j   ], XMM4
    movsd   qWORD[a_j+2*FloatSize], XMM5
{$ELSE}
begin
  a_j.x := a_j.x+b_ij.x.x*C_i.x+b_ij.y.x*C_i.y+b_ij.z.x*C_i.z ;
  a_j.y := a_j.y+b_ij.x.y*C_i.x+b_ij.y.y*C_i.y+b_ij.z.y*C_i.z ;
  a_j.z := a_j.z+b_ij.x.z*C_i.x+b_ij.y.z*C_i.y+b_ij.z.z*C_i.z ;
{$ENDIF CPUX64}

end;

class procedure T_SSE.Add(var a_j: T_Vect; const C_i: T_VectArr; const b_ij:
    T_TensorArr; const N: integer);
 {$IFDEF CPUX64}
asm
{  a_i.x := a_i.x+b_ij.x.x*C_j.x+b_ij.x.y*C_j.y+b_ij.x.z*C_j.z ;
  a_i.y := a_i.y+b_ij.y.x*C_j.x+b_ij.y.y*C_j.y+b_ij.y.z*C_j.z ;
  a_i.z := a_i.z+b_ij.z.x*C_j.x+b_ij.z.y*C_j.y+b_ij.z.z*C_j.z ;}
  // C_i: RDX,
  // b_ij:  r8   (?)
  // a_j: RCX


    .NOFRAME
    .PUSHNV RBX
    .SAVENV XMM5

    XOR RBX, RBX

    mov EBX, N


  movapd  XMM4, oWORD[ a_j              ]   // XMM4 = a0 | a1
  movlpd  XMM5, qWORD[ a_j+2*FloatSize  ]   // XMM5 = a2 | w

  @next:



    movddup   XMM7, qWORD[ C_i   ]      //  XMM7 = c0 | c0
    ADD     C_i, FloatSize                   //  C_i = C_i + 1 elem

    movapd  XMM0, oWORD[ b_ij  ]       //  XMM0 = m00 | m01
    ADD b_ij, 2*FloatSize              //  b_ij = @M+2=  @M.02
    mulpd XMM0, XMM7                 //  xmm0 = m00*c0 | m01*c0


    movlpd   XMM1, qWORD[ b_ij  ]       //  XMM1 = m02 | w
    ADD b_ij, 2*FloatSize              //  b_ij = @M+4   @M.10
    mulsd XMM1, XMM7                 //  xmm1 = m02*c0 | 0


    movddup   XMM7, qWORD[ C_i   ]    // XMM7 = c1 | c1
    ADD     C_i, FloatSize                   //  C_i = C_i + 2 elem

    addpd  XMM4, XMM0                //  a0+c0*m00 | a1 + c0*m01
    addsd  XMM5, XMM1                //  a2+c0*m02 | 0

    movapd  XMM0, oWORD[ b_ij  ]       //  XMM0 = m10 | m11
    ADD b_ij, 2*FloatSize              //  b_ij = @M+6 =@M.12
    movlpd   XMM1, qWORD[ b_ij  ]       //  XMM1 = m12 | w
    ADD b_ij, 2*FloatSize              //  b_ij = @M+8 =@M.20

    mulpd XMM0, XMM7                 // xmm0 = m10*c1 | m11*c1

    addpd  XMM4, XMM0                //  a0+c0*m00+m10*c1 | a1 + c0*m01+m11*c1


    mulsd XMM1, XMM7                 // xmm1 = m12*c1 | 0

    addsd  XMM5, XMM1                //  a2+c0*m02+m12*c1 | 0

    movddup   XMM7, qWORD[ C_i  ]     // XMM7 = c2 | c2
    ADD     C_i, 2*FloatSize                   //  C_i = C_i + 3 elem


    movapd  XMM0, oWORD[ b_ij  ]       //  XMM0 = m20 | m21
    ADD b_ij, 2*FloatSize              // b_ij = @M+10 = @M.22
    mulpd XMM0, XMM7                 // xmm0 = m20*c2 | m21*c2

    movlpd   XMM1, qWORD[ b_ij ]        //  XMM1 = m22 | w
    ADD b_ij, 2*FloatSize              // b_ij = @M+12 = @M.next ! (12= 3 Vect*4 elems)
    mulsd XMM1, XMM7                 // xmm1 = m22*c2 | 0

    addpd  XMM4, XMM0                //  a0+c0*m00+m10*c1+m20*c2 | a1 + c0*m01+m11*c1+m21*c2
    addsd  XMM5, XMM1                //  a2+c0*m02+m12*c1+m22*c2 | 0




    dec     EBX             // уменьшаем коунтер
    jnz     @next          // прыг на новую итерацию.

    movapd  oWORD[a_j              ], XMM4
    movsd   qWORD[a_j + 2*FloatSize], XMM5


{$ELSE}

var i : integer;
begin

  for i := 0 to N-1 do
  begin
    T_VectArr(@a_j)[i].x := T_VectArr(@a_j)[i].x +
       T_VectArr(@C_i)[i].x * T_TensorArr(@b_ij)[i].x.x +
       T_VectArr(@C_i)[i].y * T_TensorArr(@b_ij)[i].y.x +
       T_VectArr(@C_i)[i].z * T_TensorArr(@b_ij)[i].z.x
       ;

    T_VectArr(@a_j)[i].y := T_VectArr(@a_j)[i].y +
       T_VectArr(@C_i)[i].x * T_TensorArr(@b_ij)[i].x.y +
       T_VectArr(@C_i)[i].y * T_TensorArr(@b_ij)[i].y.y +
       T_VectArr(@C_i)[i].z * T_TensorArr(@b_ij)[i].z.y
       ;

    T_VectArr(@a_j)[i].z := T_VectArr(@a_j)[i].z +
       T_VectArr(@C_i)[i].x * T_TensorArr(@b_ij)[i].x.z +
       T_VectArr(@C_i)[i].y * T_TensorArr(@b_ij)[i].y.z +
       T_VectArr(@C_i)[i].z * T_TensorArr(@b_ij)[i].z.z
       ;

  end;


{$ENDIF CPUX64}

end;

class procedure T_SSE.Add(var a_ij: T_M4; const b_ij:
    T_M4; const c: real);
{$IFDEF CPUX64}
asm
// TEST COVERED
  .NOFRAME


  movapd  XMM4, c  // calling convention set as "Register", "C" stored in XMM2 as VALUE!!!
                   // "c" is not Pointer!!!!
  movlhps XMM4, XMM4


  movapd XMM0, oWORD[b_ij + 0*FloatSize ]
  movapd XMM1, oWORD[b_ij + 2*FloatSize ]
  movapd XMM2, oWORD[b_ij + 4*FloatSize ]
  movapd XMM3, oWORD[b_ij + 6*FloatSize ]

  mulpd  XMM0, XMM4
  mulpd  XMM1, XMM4
  mulpd  XMM2, XMM4
  mulpd  XMM3, XMM4

  addpd  XMM0, oWORD[a_ij + 0*FloatSize ]
  addpd  XMM1, oWORD[a_ij + 2*FloatSize ]
  addpd  XMM2, oWORD[a_ij + 4*FloatSize ]
  addpd  XMM3, oWORD[a_ij + 6*FloatSize ]

  movapd oWORD[a_ij + 0*FloatSize ], XMM0
  movapd oWORD[a_ij + 2*FloatSize ], XMM1
  movapd oWORD[a_ij + 4*FloatSize ], XMM2
  movapd oWORD[a_ij + 6*FloatSize ], XMM3

  movapd XMM0, oWORD[b_ij +  8*FloatSize ]
  movapd XMM1, oWORD[b_ij + 10*FloatSize ]
  movapd XMM2, oWORD[b_ij + 12*FloatSize ]
  movapd XMM3, oWORD[b_ij + 14*FloatSize ]

  mulpd  XMM0, XMM4
  mulpd  XMM1, XMM4
  mulpd  XMM2, XMM4
  mulpd  XMM3, XMM4

  addpd  XMM0, oWORD[a_ij + 8*FloatSize ]
  addpd  XMM1, oWORD[a_ij + 10*FloatSize ]
  addpd  XMM2, oWORD[a_ij + 12*FloatSize ]
  addpd  XMM3, oWORD[a_ij + 14*FloatSize ]

  movapd oWORD[a_ij +  8*FloatSize ], XMM0
  movapd oWORD[a_ij + 10*FloatSize ], XMM1
  movapd oWORD[a_ij + 12*FloatSize ], XMM2
  movapd oWORD[a_ij + 14*FloatSize ], XMM3

{$ELSE} {X32 MODE}  begin
  a_ij[ 0, 0] := a_ij[ 0, 0] + b_ij[ 0, 0] *c;
  a_ij[ 0, 1] := a_ij[ 0, 1] + b_ij[ 0, 1] *c;
  a_ij[ 0, 2] := a_ij[ 0, 2] + b_ij[ 0, 2] *c;
  a_ij[ 0, 3] := a_ij[ 0, 3] + b_ij[ 0, 3] *c;

  a_ij[ 1, 0] := a_ij[ 1, 0] + b_ij[ 1, 0] *c;
  a_ij[ 1, 1] := a_ij[ 1, 1] + b_ij[ 1, 1] *c;
  a_ij[ 1, 2] := a_ij[ 1, 2] + b_ij[ 1, 2] *c;
  a_ij[ 1, 3] := a_ij[ 1, 3] + b_ij[ 1, 3] *c;

  a_ij[ 2, 0] := a_ij[ 2, 0] + b_ij[ 2, 0] *c;
  a_ij[ 2, 1] := a_ij[ 2, 1] + b_ij[ 2, 1] *c;
  a_ij[ 2, 2] := a_ij[ 2, 2] + b_ij[ 2, 2] *c;
  a_ij[ 2, 3] := a_ij[ 2, 3] + b_ij[ 2, 3] *c;

  a_ij[ 3, 0] := a_ij[ 3, 0] + b_ij[ 3, 0] *c;
  a_ij[ 3, 1] := a_ij[ 3, 1] + b_ij[ 3, 1] *c;
  a_ij[ 3, 2] := a_ij[ 3, 2] + b_ij[ 3, 2] *c;
  a_ij[ 3, 3] := a_ij[ 3, 3] + b_ij[ 3, 3] *c;
{$ENDIF}
end;

class procedure T_SSE.Add(var a_ij: T_M4; const DC_i, DC_j:
    T_Vect4);
{$IFDEF CPUX64}
asm
// TEST COVERED
  .NOFRAME
  .SAVENV XMM4
  .SAVENV XMM5


  movapd XMM0, oWORD[ DC_i + 0*FloatSize ]
  movapd XMM1, oWORD[ DC_i + 2*FloatSize ]

  movapd XMM2, oWORD[ DC_j + 0*FloatSize ]
  movapd XMM3, oWORD[ DC_j + 2*FloatSize ]

  ///// a_ij[-1]
  movapd  XMM4, XMM0
  movlhps XMM4, XMM4 // XMM4 = DC_i[-1] |  DC_i[-1]

  movapd  XMM5, XMM4 // XMM5 = DC_i[-1] |  DC_i[-1]

  mulpd   XMM4, XMM2 // XMM4 = DC_i[-1] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = DC_i[-1] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[a_ij +  0*FloatSize ]
  addpd XMM5, oWORD[a_ij +  2*FloatSize ]

  movapd oWORD[a_ij +  0*FloatSize ], XMM4
  movapd oWORD[a_ij +  2*FloatSize ], XMM5

  ///// Next stringL a_ij[0]
  ADD a_ij, 4*FloatSize

  movapd  XMM4, XMM0
  movhlps XMM4, XMM4 // XMM4 = DC_i[ 0] |  DC_i[ 0]

  movapd  XMM5, XMM4 // XMM5 = DC_i[ 0] |  DC_i[ 0]
  mulpd   XMM4, XMM2 // XMM4 = DC_i[ 0] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = DC_i[ 0] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[a_ij +  0*FloatSize ]
  addpd XMM5, oWORD[a_ij +  2*FloatSize ]

  movapd oWORD[a_ij +  0*FloatSize ], XMM4
  movapd oWORD[a_ij +  2*FloatSize ], XMM5

  ///// Next string  a_ij[1]
  ADD a_ij, 4*FloatSize

  movapd  XMM4, XMM1
  movlhps XMM4, XMM4 // XMM4 = DC_i[ 1] |  DC_i[ 1]

  movapd  XMM5, XMM4 // XMM5 = DC_i[ 1] |  DC_i[ 1]
  mulpd   XMM4, XMM2 // XMM4 = DC_i[ 1] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = DC_i[ 1] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[a_ij +  0*FloatSize ]
  addpd XMM5, oWORD[a_ij +  2*FloatSize ]

  movapd oWORD[a_ij +  0*FloatSize ], XMM4
  movapd oWORD[a_ij +  2*FloatSize ], XMM5

  ///// Next string   a_ij[2]
  ADD a_ij, 4*FloatSize

  movapd  XMM4, XMM1
  movhlps XMM4, XMM4 // XMM4 = DC_i[ 2] |  DC_i[ 2]

  movapd  XMM5, XMM4 // XMM5 = DC_i[ 2] |  DC_i[ 2]
  mulpd   XMM4, XMM2 // XMM4 = DC_i[ 2] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = DC_i[ 2] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[a_ij +  0*FloatSize ]
  addpd XMM5, oWORD[a_ij +  2*FloatSize ]

  movapd oWORD[a_ij +  0*FloatSize ], XMM4
  movapd oWORD[a_ij +  2*FloatSize ], XMM5


{$ELSE} {X32 MODE}begin
  a_ij[ 0, 0] := a_ij[ 0, 0] + DC_j[ 0]*DC_i[ 0];
  a_ij[ 0, 1] := a_ij[ 0, 1] + DC_j[ 1]*DC_i[ 0];
  a_ij[ 0, 2] := a_ij[ 0, 2] + DC_j[ 2]*DC_i[ 0];
  a_ij[ 0, 3] := a_ij[ 0, 3] + DC_j[ 3]*DC_i[ 0];

  a_ij[ 1, 0] := a_ij[ 1, 0] + DC_j[ 0]*DC_i[ 1];
  a_ij[ 1, 1] := a_ij[ 1, 1] + DC_j[ 1]*DC_i[ 1];
  a_ij[ 1, 2] := a_ij[ 1, 2] + DC_j[ 2]*DC_i[ 1];
  a_ij[ 1, 3] := a_ij[ 1, 3] + DC_j[ 3]*DC_i[ 1];

  a_ij[ 2, 0] := a_ij[ 2, 0] + DC_j[ 0]*DC_i[ 2];
  a_ij[ 2, 1] := a_ij[ 2, 1] + DC_j[ 1]*DC_i[ 2];
  a_ij[ 2, 2] := a_ij[ 2, 2] + DC_j[ 2]*DC_i[ 2];
  a_ij[ 2, 3] := a_ij[ 2, 3] + DC_j[ 3]*DC_i[ 2];

  a_ij[ 3, 0] := a_ij[ 3, 0] + DC_j[ 0]*DC_i[ 3];
  a_ij[ 3, 1] := a_ij[ 3, 1] + DC_j[ 1]*DC_i[ 3];
  a_ij[ 3, 2] := a_ij[ 3, 2] + DC_j[ 2]*DC_i[ 3];
  a_ij[ 3, 3] := a_ij[ 3, 3] + DC_j[ 3]*DC_i[ 3];
{$ENDIF}
end;

class procedure T_SSE.Add(var DC_i: T_Vect4; const a_ij:
    T_M4; const DC_j: T_Vect4);
{$IFDEF CPUX64}
asm
// TEST COVERED
  .NOFRAME

  .SAVENV XMM5

  movapd XMM0, oWORD[  DC_j + 0*FloatSize ]
  movapd XMM1, oWORD[  DC_j + 2*FloatSize ]

  movapd XMM2, oWORD[  a_ij + 0*FloatSize ]
  movapd XMM3, oWORD[  a_ij + 2*FloatSize ]

  movapd XMM4, oWORD[  a_ij + 4*FloatSize ]
  movapd XMM5, oWORD[  a_ij + 6*FloatSize ]

  mulpd XMM2, XMM0
  mulpd XMM3, XMM1

  mulpd XMM4, XMM0
  mulpd XMM5, XMM1

  haddpd XMM2,XMM3
  haddpd XMM4,XMM5

  haddpd XMM2, XMM4


  addpd XMM2, oWORD [ DC_i + 0*FloatSize]

  movapd oWORD [ DC_i + 0*FloatSize], XMM2




  movapd XMM2, oWORD[  a_ij +  8*FloatSize ]
  movapd XMM3, oWORD[  a_ij + 10*FloatSize ]

  movapd XMM4, oWORD[  a_ij + 12*FloatSize ]
  movapd XMM5, oWORD[  a_ij + 14*FloatSize ]

  mulpd XMM2, XMM0
  mulpd XMM3, XMM1

  mulpd XMM4, XMM0
  mulpd XMM5, XMM1

  haddpd XMM2,XMM3
  haddpd XMM4,XMM5

  haddpd XMM2, XMM4


  addpd XMM2, oWORD [ DC_i + 2*FloatSize]

  movapd oWORD [ DC_i + 2*FloatSize], XMM2


{$ELSE} {X32 MODE}
  begin
    DC_i[-1] :=  DC_i[-1] +
      a_ij[-1,-1] *DC_j[-1]+
      a_ij[-1, 0] *DC_j[ 0]+
      a_ij[-1, 1] *DC_j[ 1]+
      a_ij[-1, 2] *DC_j[ 2];


    DC_i[ 0] := DC_i[ 0] +
      a_ij[ 0,-1] *DC_j[-1]+
      a_ij[ 0, 0] *DC_j[ 0]+
      a_ij[ 0, 1] *DC_j[ 1]+
      a_ij[ 0, 2] *DC_j[ 2];


    DC_i[ 1] := DC_i[ 1] + a_ij[ 1,-1] *DC_j[-1]+a_ij[ 1, 0] *DC_j[ 0]+a_ij[ 1, 1] *DC_j[ 1]+a_ij[ 1, 2] *DC_j[ 2];
    DC_i[ 2] := DC_i[ 2] + a_ij[ 2,-1] *DC_j[-1]+a_ij[ 2, 0] *DC_j[ 0]+a_ij[ 2, 1] *DC_j[ 1]+a_ij[ 2, 2] *DC_j[ 2];

{$ENDIF}
end;
class procedure T_SSE.Add(var a: T_Vect4; const b: T_Vect4);
begin
  a[0] := a[0]+b[0];
  a[1] := a[1]+b[1];
  a[2] := a[2]+b[2];
  a[3] := a[3]+b[3];
end;

class procedure T_SSE.Add(var a: T_Vect4; const b: T_Vect4; const
    koeff: real);
begin
  a[ 0] := a[ 0]+b[ 0]*koeff;
  a[ 1] := a[ 1]+b[ 1]*koeff;
  a[ 2] := a[ 2]+b[ 2]*koeff;
  a[ 3] := a[ 3]+b[ 3]*koeff;
end;

class procedure T_SSE.Add(var a: T_Vect4; const b: T_Vect; const koeff:
    real);
begin
  a[ 0] := a[ 0]+b.X*koeff;
  a[ 1] := a[ 1]+b.Y*koeff;
  a[ 2] := a[ 2]+b.Z*koeff;
end;

class procedure T_SSE.Add(var a: T_Tens; const b: T_Tens);
begin
  a.x.x := a.x.x + b.x.x;
  a.y.x := a.y.x + b.y.x;
  a.z.x := a.z.x + b.z.x;

  a.x.y := a.x.y + b.x.y;
  a.y.y := a.y.y + b.y.y;
  a.z.y := a.z.y + b.z.y;

  a.x.z := a.x.z + b.x.z;
  a.y.z := a.y.z + b.y.z;
  a.z.z := a.z.z + b.z.z;
end;

class procedure T_SSE.Add(var a: T_Tens; const b: T_Tens; const koeff: real);
begin
  a.x.x := a.x.x + b.x.x*koeff;
  a.x.y := a.x.y + b.x.y*koeff;
  a.x.z := a.x.z + b.x.z*koeff;

  a.y.x := a.y.x + b.y.x*koeff;
  a.y.y := a.y.y + b.y.y*koeff;
  a.y.z := a.y.z + b.y.z*koeff;

  a.z.x := a.z.x + b.z.x*koeff;
  a.z.y := a.z.y + b.z.y*koeff;
  a.z.z := a.z.z + b.z.z*koeff;
end;

class procedure T_SSE.Add(var a: T_Tens; const b: real);
begin
  a.x.x := a.x.x + b;
  a.y.y := a.y.y + b;
  a.z.z := a.z.z + b;
end;

class procedure T_SSE.Add(var B_ab: T_Tens; const _DR_alfa, _DR_beta: T_Vect);
{$IFDEF CPUX64}
asm
// TEST COVERED
  .NOFRAME
  .SAVENV XMM4
  .SAVENV XMM5
  .SAVENV XMM6

  movapd XMM0, oWORD[ _DR_alfa + 0*FloatSize ]
  movapd XMM1, oWORD[ _DR_alfa + 2*FloatSize ]

  movapd XMM2, oWORD[ _DR_beta + 0*FloatSize ]
  movapd XMM3, oWORD[ _DR_beta + 2*FloatSize ]

  ///// a_ij[-1]
  movapd  XMM4, XMM0
  movlhps XMM4, XMM4 // XMM4 = _DR_alfa[-1] |  _DR_alfa[-1]

  movapd  XMM5, XMM4 // XMM5 = _DR_alfa[-1] |  _DR_alfa[-1]

  mulpd   XMM4, XMM2 // XMM4 = _DR_alfa[-1] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = _DR_alfa[-1] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[B_ab +  0*FloatSize ]
  addpd XMM5, oWORD[B_ab +  2*FloatSize ]

  movapd oWORD[B_ab +  0*FloatSize ], XMM4
  movlpd qWORD[B_ab +  2*FloatSize ], XMM5

  ///// Next stringL a_ij[0]
  ADD B_ab, 4*FloatSize

  movapd  XMM4, XMM0
  movhlps XMM4, XMM4 // XMM4 = _DR_alfa[ 0] |  _DR_alfa[ 0]

  movapd  XMM5, XMM4 // XMM5 = _DR_alfa[ 0] |  _DR_alfa[ 0]
  mulpd   XMM4, XMM2 // XMM4 = _DR_alfa[ 0] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = _DR_alfa[ 0] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[B_ab +  0*FloatSize ]
  addpd XMM5, oWORD[B_ab +  2*FloatSize ]

  movapd oWORD[B_ab +  0*FloatSize ], XMM4
  movlpd qWORD[B_ab +  2*FloatSize ], XMM5

  ///// Next string  a_ij[1]
  ADD B_ab, 4*FloatSize

  movapd  XMM4, XMM1
  movlhps XMM4, XMM4 // XMM4 = _DR_alfa[ 1] |  _DR_alfa[ 1]

  movapd  XMM5, XMM4 // XMM5 = _DR_alfa[ 1] |  DC_i[ 1]
  mulpd   XMM4, XMM2 // XMM4 = _DR_alfa[ 1] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = _DR_alfa[ 1] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[B_ab +  0*FloatSize ]
  addpd XMM5, oWORD[B_ab +  2*FloatSize ]

  movapd oWORD[B_ab +  0*FloatSize ], XMM4
  movlpd qWORD[B_ab +  2*FloatSize ], XMM5


{$ELSE} {X32 MODE}begin
    // B[alfa, beta ] = _DR[alfa]*DC_beta,  dR[-1] = 1.0;

    B_ab.x.x := B_ab.x.x + _DR_alfa.x * _DR_beta.x;
    B_ab.x.y := B_ab.x.y + _DR_alfa.x * _DR_beta.y;
    B_ab.x.z := B_ab.x.z + _DR_alfa.x * _DR_beta.z;


    B_ab.y.x := B_ab.y.x + _DR_alfa.y * _DR_beta.x;
    B_ab.y.y := B_ab.y.y + _DR_alfa.y * _DR_beta.y;
    B_ab.y.z := B_ab.y.z + _DR_alfa.y * _DR_beta.z;


    B_ab.z.x := B_ab.z.x + _DR_alfa.z * _DR_beta.x;
    B_ab.z.y := B_ab.z.y + _DR_alfa.z * _DR_beta.y;
    B_ab.z.z := B_ab.z.z + _DR_alfa.z * _DR_beta.z;
{$ENDIF}
end;

class procedure T_SSE.Add(var B_ab: T_Tens; const _DR_alfa, _DR_beta: T_Vect;
    const alfa: real);
{$IFDEF CPUX64}
asm
// TEST COVERED

  .NOFRAME
  .SAVENV XMM4
  .SAVENV XMM5
  .SAVENV XMM6

//  movlpd  XMM4, qWORD[alfa]
  movapd  XMM4, alfa // Alfa stored in XMM3! becouse of "register" calling convention
  movlhps XMM4, XMM4
  movapd XMM0, oWORD[ _DR_alfa + 0*FloatSize ]
  movapd XMM1, oWORD[ _DR_alfa + 2*FloatSize ]

  mulpd XMM0, XMM4 //  _DR_alfa* Coeff
  mulsd XMM1, XMM4 //

  movapd XMM2, oWORD[ _DR_beta + 0*FloatSize ]
  movapd XMM3, oWORD[ _DR_beta + 2*FloatSize ]

  ///// a_ij[-1]
  movapd  XMM4, XMM0
  movlhps XMM4, XMM4 // XMM4 = _DR_alfa[-1] |  _DR_alfa[-1]

  movapd  XMM5, XMM4 // XMM5 = _DR_alfa[-1] |  _DR_alfa[-1]

  mulpd   XMM4, XMM2 // XMM4 = _DR_alfa[-1] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = _DR_alfa[-1] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[B_ab +  0*FloatSize ]
  addpd XMM5, oWORD[B_ab +  2*FloatSize ]

  movapd oWORD[B_ab +  0*FloatSize ], XMM4
  movlpd qWORD[B_ab +  2*FloatSize ], XMM5

  ///// Next stringL a_ij[0]
  ADD B_ab, 4*FloatSize

  movapd  XMM4, XMM0
  movhlps XMM4, XMM4 // XMM4 = _DR_alfa[ 0] |  _DR_alfa[ 0]

  movapd  XMM5, XMM4 // XMM5 = _DR_alfa[ 0] |  _DR_alfa[ 0]
  mulpd   XMM4, XMM2 // XMM4 = _DR_alfa[ 0] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = _DR_alfa[ 0] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[B_ab +  0*FloatSize ]
  addpd XMM5, oWORD[B_ab +  2*FloatSize ]

  movapd oWORD[B_ab +  0*FloatSize ], XMM4
  movlpd qWORD[B_ab +  2*FloatSize ], XMM5

  ///// Next string  a_ij[1]
  ADD B_ab, 4*FloatSize

  movapd  XMM4, XMM1
  movlhps XMM4, XMM4 // XMM4 = _DR_alfa[ 1] |  _DR_alfa[ 1]

  movapd  XMM5, XMM4 // XMM5 = _DR_alfa[ 1] |  DC_i[ 1]
  mulpd   XMM4, XMM2 // XMM4 = _DR_alfa[ 1] * (DC_j[-1] | DC_j[ 0])
  mulpd   XMM5, XMM3 // XMM5 = _DR_alfa[ 1] * (DC_j[ 1] | DC_j[ 2])

  addpd XMM4, oWORD[B_ab +  0*FloatSize ]
  addpd XMM5, oWORD[B_ab +  2*FloatSize ]

  movapd oWORD[B_ab +  0*FloatSize ], XMM4
  movlpd qWORD[B_ab +  2*FloatSize ], XMM5


{$ELSE} {X32 MODE}begin
    // B[alfa, beta ] = _DR[alfa]*DC_beta,  dR[-1] = 1.0;

    B_ab.x.x := B_ab.x.x + _DR_alfa.x * _DR_beta.x*alfa;
    B_ab.x.y := B_ab.x.y + _DR_alfa.x * _DR_beta.y*alfa;
    B_ab.x.z := B_ab.x.z + _DR_alfa.x * _DR_beta.z*alfa;


    B_ab.y.x := B_ab.y.x + _DR_alfa.y * _DR_beta.x*alfa;
    B_ab.y.y := B_ab.y.y + _DR_alfa.y * _DR_beta.y*alfa;
    B_ab.y.z := B_ab.y.z + _DR_alfa.y * _DR_beta.z*alfa;


    B_ab.z.x := B_ab.z.x + _DR_alfa.z * _DR_beta.x*alfa;
    B_ab.z.y := B_ab.z.y + _DR_alfa.z * _DR_beta.y*alfa;
    B_ab.z.z := B_ab.z.z + _DR_alfa.z * _DR_beta.z*alfa;
{$ENDIF}
end;

class procedure T_SSE.Add(var a: T_Vect; const C: real; const b: T_Vect);
begin
  a.x := a.x+b.x*C;
  a.y := a.y+b.y*C;
  a.z := a.z+b.z*C;
end;

class procedure T_SSE.Add(var a: T_Vect; const C: real; const b: T_Vect;
    const d: real);
begin
  Add(a, B, C*d);
end;

class procedure T_SSE.Add(var a_i: T_Vect; const b_ij: T_Tens; const C_j:
    T_Vect4);
{$IFDEF CPUX64}
asm
// TEST COVERED
    .NOFRAME
    .SAVENV XMM5
    .SAVENV XMM6
    .SAVENV XMM7


    movapd  XMM6, oWORD[C_j]         // XMM6 = c-1 | c0
    movapd  XMM7, oWORD[C_j+16]      // XMM7 = c1  | c2

    movhlps XMM6, XMM6
    movlhps XMM6, XMM7
    movhlps XMM7, XMM7

    movupd  XMM4, oWORD[a_i]         // XMM4 = a0 | a1
    movsd  XMM5, qWORD[a_i+16]      // XMM5 = a2 | w

    movupd  XMM0, oWORD[ b_ij+0  ]  //  XMM0 = m00 | m01
    movsd  XMM1, qWORD[ b_ij+16 ]  //  XMM1 = m02 | w

    movupd  XMM2, oWORD[ b_ij+32 ]  //  XMM2 = m10 | m11
    movsd  XMM3, qWORD[ b_ij+48 ]  //  XMM3 = m12 | w


    mulpd XMM0, XMM6            // xmm0 = m00*c0 | m01*c1
    mulsd XMM1, XMM7            // xmm1 = m02*c2 | 0
    haddpd XMM0, XMM1           // xmm0 = m00*c0 + m01*c1 | m02*c2

    mulpd XMM2, XMM6            // xmm2 = m10*c0 | m11*c1
    mulsd XMM3, XMM7            // xmm3 = m22*c2 | 0
    haddpd XMM2, XMM3           // xmm2 = m10*c0 + m11*c1 | m12*c2


    haddpd XMM0, XMM2           // xmm0 = m00*c0 + m01*c1 + m02*c2 | m10*c0 + m11*c1 + m12*c2
    addpd  XMM4, XMM0

    movupd  XMM2, oWORD[ b_ij+64 ]  //  XMM2 = m20 | m21
    movsd  XMM3, qWORD[ b_ij+80 ]  //  XMM3 = m22 | w

    mulpd XMM2, XMM6            // XMM2 = m20*c0 | m21*c1
    mulsd XMM3, XMM7            // XMM3 = m22*c2 | 0
    movlhps XMM5, XMM3          // XMM5 = a2, m22*c2

    haddpd XMM5, XMM2           // xmm3 = m22*c2+a2 | m20*c0 + m21*c1
    haddpd XMM5, XMM5

    movupd  oWORD[a_i], XMM4
    movsd  qWORD[a_i+16], XMM5


{$ELSE}
begin
  a_i.x := a_i.x+ b_ij.x.x*C_j[0]+b_ij.x.y*C_j[1]+b_ij.x.z*C_j[2] ;
  a_i.y := a_i.y+ b_ij.y.x*C_j[0]+b_ij.y.y*C_j[1]+b_ij.y.z*C_j[2] ;
  a_i.z := a_i.z+ b_ij.z.x*C_j[0]+b_ij.z.y*C_j[1]+b_ij.z.z*C_j[2] ;
{$ENDIF}
end;

class procedure T_SSE.Add(var a: T_Vect; const b: T_Vect);
{$IFDEF CPUX64}
asm
  .NOFRAME

  movapd XMM0, oWORD[ b ]
  movapd XMM1, oWORD[ b + 2*FloatSize ]

  addpd XMM0, oWORD[ a]
  addsd XMM1, qWORD[ a+ 2*FloatSize]

  movapd oWORD[a             ], XMM0
  movlpd qWORD[a+ 2*FloatSize], XMM1


{$ELSE}
begin
  a.x := a.x+b.x;
  a.y := a.y+b.y;
  a.z := a.z+b.z;
{$ENDIF}
end;

class procedure T_SSE.Add(var a: T_Vect; const b: T_Vect; const C: real);
{$IFDEF CPUX64}
asm
  .NOFRAME

  movapd  XMM2, C
  movlhps XMM2, XMM2
  movapd XMM0, oWORD[ b ]
  movapd XMM1, oWORD[ b + 2*FloatSize ]
  mulpd XMM0, XMM2
  mulsd XMM1, XMM2

  addpd XMM0, oWORD[ a]
  addsd XMM1, qWORD[ a+ 2*FloatSize]

  movapd oWORD[a             ], XMM0
  movlpd qWORD[a+ 2*FloatSize], XMM1


{$ELSE}
begin
  a.x := a.x+b.x*C;
  a.y := a.y+b.y*C;
  a.z := a.z+b.z*C;
{$ENDIF}
end;

class procedure T_SSE.Add(var a: T_Vect; const b : T_Vect; const c, V: real);
{$IFDEF CPUX64}
asm
  .NOFRAME

  movapd   XMM0, c
  movapd   XMM1, V

  mulsd   XMM0, XMM1
  movlhps XMM0, XMM0

  movapd XMM1, oWORD[b + 0*FloatSize]
  movapd XMM2, oWORD[b + 2*FloatSize]

  mulpd XMM1, XMM0
  mulpd XMM2, XMM0

  addpd XMM1, oWORD[ a + 0*FloatSize ]
  addsd XMM2, qWORD[ a + 2*FloatSize ]

  movapd oWORD[ a + 0*FloatSize ], XMM1
  movlpd qWORD[ a + 2*FloatSize ], XMM2


{$ELSE}
begin
  a.x := a.x + b.x*c*V;
  a.y := a.y + b.y*c*V;
  a.z := a.z + b.z*c*V;
{$ENDIF}
end;

class procedure T_SSE.Add(var a_i: T_Vect; const C_j: T_Vect; const b_ij:
    T_Tens; const K: real);
{$IFDEF CPUX64}
asm
  .NOFRAME

  movapd   XMM2, K
  movapd   XMM0, oWORD[ C_j + 0*FloatSize ]
  movapd   XMM1, oWORD[ C_j + 2*FloatSize ]


  movlhps XMM2, XMM2
  mulpd   XMM0, XMM2
  mulpd   XMM1, XMM2


  movapd  XMM2, oWORD[ b_ij + 0*FloatSize ]
  movapd  XMM3, oWORD[ b_ij + 2*FloatSize ]
  mulpd   XMM2, XMM0    // xmm2 = b.xx*c.x | b.xy*c.y
  mulsd   XMM3, XMM1    // xmm3 = b.xz*c.z | trash
  addsd   XMM2, XMM3    // xmm2 = b.xx*c.x + b.xz*c.z| b.xy*c.y


  movapd  XMM3, oWORD[ b_ij + 4*FloatSize ]
  movapd  XMM4, oWORD[ b_ij + 6*FloatSize ]
  mulpd   XMM3, XMM0    // xmm3 = b.yx*c.x | b.yy*c.y
  mulsd   XMM4, XMM1    // xmm4 = b.yz*c.z | trash
  addsd   XMM3, XMM4    // xmm3 = b.yx*c.x + b.yz*c.z| b.yy*c.y
  haddpd  XMM2, XMM3

  addpd   XMM2, oWORD[ a_i + 0*FloatSize]
  movapd oWORD[ a_i + 0*FloatSize ], XMM2

  movapd  XMM2, oWORD[ b_ij +  8*FloatSize ]
  movapd  XMM3, oWORD[ b_ij + 10*FloatSize ]
  mulpd   XMM2, XMM0    // xmm2 = b.zx*c.x | b.zy*c.y
  mulsd   XMM3, XMM1    // xmm3 = b.zz*c.z | trash
  movhpd  XMM3, qWORD [a_i + 2*FloatSize] //// xmm3 = b.zz*c.z | a_i.z

  haddpd XMM2, XMM3  // XMM2 =  b.zx*c.x + b.zy*c.y | b.zz*c.z + a_i.z
  haddpd XMM2, XMM2  // XMM2 =  b.zx*c.x + b.zy*c.y + b.zz*c.z + a_i.z  | twice
  movlpd  qWORD [a_i + 2*FloatSize], XMM2

{$ELSE}
begin
  a_i.x := a_i.x+K*(b_ij.x.x*C_j.x+b_ij.x.y*C_j.y+b_ij.x.z*C_j.z) ;
  a_i.y := a_i.y+K*(b_ij.y.x*C_j.x+b_ij.y.y*C_j.y+b_ij.y.z*C_j.z) ;
  a_i.z := a_i.z+K*(b_ij.z.x*C_j.x+b_ij.z.y*C_j.y+b_ij.z.z*C_j.z) ;
{$ENDIF}
end;

class procedure T_SSE.Add(const a_j: T_VectArr; const C_i: T_VectArr; const
    b_ij: T_TensorArr; const N: integer);
 {$IFDEF CPUX64}
asm
{  a_i.x := a_i.x+b_ij.x.x*C_j.x+b_ij.x.y*C_j.y+b_ij.x.z*C_j.z ;
  a_i.y := a_i.y+b_ij.y.x*C_j.x+b_ij.y.y*C_j.y+b_ij.y.z*C_j.z ;
  a_i.z := a_i.z+b_ij.z.x*C_j.x+b_ij.z.y*C_j.y+b_ij.z.z*C_j.z ;}
  // C_i: RDX,
  // b_ij:  r8   (?)
  // a_j: RCX


    .NOFRAME
    .PUSHNV RBX
    .SAVENV XMM5

    XOR RBX, RBX

    mov EBX, N



  @next:

    movupd  XMM4, oWORD[ a_j              ]   // XMM4 = a0 | a1
    movlpd  XMM5, qWORD[ a_j+2*FloatSize  ]   // XMM5 = a2 | w

    movddup   XMM7, qWORD[ C_i   ]      //  XMM7 = c0 | c0
    Add     C_i, FloatSize                   //  C_i = C_i + 1 elem

    movupd  XMM0, oWORD[ b_ij  ]       //  XMM0 = m00 | m01
    Add b_ij, 2*FloatSize              //  b_ij = @M+2=  @M.02
    mulpd XMM0, XMM7                 //  xmm0 = m00*c0 | m01*c0


    movlpd   XMM1, qWORD[ b_ij  ]       //  XMM1 = m02 | w
    Add b_ij, 2*FloatSize              //  b_ij = @M+4   @M.10
    mulsd XMM1, XMM7                 //  xmm1 = m02*c0 | 0


    movddup   XMM7, qWORD[ C_i   ]    // XMM7 = c1 | c1
    Add     C_i, FloatSize                   //  C_i = C_i + 2 elem

    addpd  XMM4, XMM0                //  a0+c0*m00 | a1 + c0*m01
    addsd  XMM5, XMM1                //  a2+c0*m02 | 0

    movupd  XMM0, oWORD[ b_ij  ]       //  XMM0 = m10 | m11
    Add b_ij, 2*FloatSize              //  b_ij = @M+6 =@M.12
    movlpd   XMM1, qWORD[ b_ij  ]       //  XMM1 = m12 | w
    Add b_ij, 2*FloatSize              //  b_ij = @M+8 =@M.20

    mulpd XMM0, XMM7                 // xmm0 = m10*c1 | m11*c1

    addpd  XMM4, XMM0                //  a0+c0*m00+m10*c1 | a1 + c0*m01+m11*c1


    mulsd XMM1, XMM7                 // xmm1 = m12*c1 | 0

    addsd  XMM5, XMM1                //  a2+c0*m02+m12*c1 | 0

    movddup   XMM7, qWORD[ C_i  ]     // XMM7 = c2 | c2
    Add     C_i, 2*FloatSize                   //  C_i = C_i + 3 elem


    movupd  XMM0, oWORD[ b_ij  ]       //  XMM0 = m20 | m21
    Add b_ij, 2*FloatSize              // b_ij = @M+10 = @M.22
    mulpd XMM0, XMM7                 // xmm0 = m20*c2 | m21*c2

    movlpd   XMM1, qWORD[ b_ij ]        //  XMM1 = m22 | w
    Add b_ij, 2*FloatSize              // b_ij = @M+12 = @M.next ! (12= 3 Vect*4 elems)
    mulsd XMM1, XMM7                 // xmm1 = m22*c2 | 0

    addpd  XMM4, XMM0                //  a0+c0*m00+m10*c1+m20*c2 | a1 + c0*m01+m11*c1+m21*c2
    addsd  XMM5, XMM1                //  a2+c0*m02+m12*c1+m22*c2 | 0

    movupd  oWORD[a_j              ], XMM4
    movsd   qWORD[a_j + 2*FloatSize], XMM5

    Add     a_j, 4*FloatSize


    dec     EBX             // уменьшаем коунтер
    jg     @next          // прыг на новую итерацию.

    Add     a_j, 4*FloatSize

{$ELSE}

var i : integer;
begin

  for i := 0 to N-1 do
  begin
    T_VectArr(@a_j)[i].x := T_VectArr(@a_j)[i].x +
       T_VectArr(@C_i)[i].x * T_TensorArr(@b_ij)[i].x.x +
       T_VectArr(@C_i)[i].y * T_TensorArr(@b_ij)[i].y.x +
       T_VectArr(@C_i)[i].z * T_TensorArr(@b_ij)[i].z.x
       ;

    T_VectArr(@a_j)[i].y := T_VectArr(@a_j)[i].y +
       T_VectArr(@C_i)[i].x * T_TensorArr(@b_ij)[i].x.y +
       T_VectArr(@C_i)[i].y * T_TensorArr(@b_ij)[i].y.y +
       T_VectArr(@C_i)[i].z * T_TensorArr(@b_ij)[i].z.y
       ;

    T_VectArr(@a_j)[i].z := T_VectArr(@a_j)[i].z +
       T_VectArr(@C_i)[i].x * T_TensorArr(@b_ij)[i].x.z +
       T_VectArr(@C_i)[i].y * T_TensorArr(@b_ij)[i].y.z +
       T_VectArr(@C_i)[i].z * T_TensorArr(@b_ij)[i].z.z
       ;

  end;


{$ENDIF CPUX64}

end;

class procedure T_SSE.Invert(var S: T_Tens);
{$IFDEF CPUX64}
  const FloatSize = 8;
        One: double = 1.0;
        Zero: double = 0.0;
        eps: double = 1.0e-10;
        sign: Uint64 = $7FFFFFFFFFFFFFFF;
  asm
    .SAVENV XMM5   // S.z.z
    .SAVENV XMM6   // Res.x.x | Res.x.y
    .SAVENV XMM7   // Res.x.z | w
    .SAVENV XMM8   // Res.y.x | Res.y.y
    .SAVENV XMM9   // Res.y.z | w
    .SAVENV XMM10  // Res.z.x | Res.z.y
    .SAVENV XMM11  // Res.z.z | w
    .SAVENV XMM12  // tmp1
    .SAVENV XMM13  // tmp2        }

    prefetch [RCX]
    prefetch [RCX+32]
    prefetch [RCX+64]


    movupd XMM0, [RCX              ]   //   S.x.x | S.x.y
    movupd XMM1, [RCX+  2*FloatSize]   //   S.x.z | w
    movupd XMM2, [RCX+  4*FloatSize]   //   S.y.x | S.y.y
    movupd XMM3, [RCX+  6*FloatSize]   //   S.y.z | w
    movupd XMM4, [RCX+  8*FloatSize]   //   S.z.x | S.z.y
    movupd XMM5, [RCX+ 10*FloatSize]   //   S.z.z | w

  {   Res.x.x :=    ( + S.y.y*S.z.z - S.y.z*S.z.y  );   XMM6
      Res.x.y :=    ( + S.x.z*S.z.y - S.x.y*S.z.z  );          }

    movhlps XMM6, XMM2  // XMM6.L = S.y.y
    movlhps XMM6, XMM3  // XMM6.h = S.y.z
    blendpd XMM5, XMM4, 2  // XMM5= S.z.z | S.z.y

    mulpd XMM6, XMM5   // XMM6 = S.y.y*S.z.z |  S.y.z*S.z.y


    movhlps XMM12, XMM4  // XMM12.L = S.z.y
    movlhps XMM12, XMM5  // XMM12.h = S.z.z
    blendpd XMM1,  XMM0, 2 // XMM1 = S.x.z | S.x.y

    mulpd XMM12, XMM1   // XMM13 = S.x.z*S.z.y | S.x.y* S.z.z

    hsubpd XMM6, XMM12   // XMM6 = Res.x.x | Res.x.y     }

  {  Res.y.x :=    ( + S.z.x*S.y.z - S.y.x*S.z.z  );
     Res.y.y :=    ( + S.x.x*S.z.z - S.x.z*S.z.x  );   }

    movsd   XMM8, XMM4  // XMM8.L = S.z.x
    movlhps XMM8, XMM2  // XMM8.h = S.y.x
  //  movsd   XMM12, XMM3  // XMM12 = S.y.z | S.z.z (ALREDY STORED!!!)
    movlhps XMM3 , XMM5  // XMM3 = S.y.z | S.z.z

    mulpd   XMM8, XMM3   // XMM8 = S.z.x*S.y.z |  S.y.x*S.z.z

    movsd   XMM13, XMM0  // XMM12.L = S.x.x
    movlhps XMM13, XMM1  // XMM13 = S.x.x | S.x.z

    movlhps XMM5, XMM4   // XMM5 = S.z.z  | S.z.x

    mulpd XMM13, XMM5    // XMM13 = S.x.x*S.z.z | S.x.z* S.z.x

    hsubpd XMM8, XMM13   // XMM8 = Res.y.x | Res.y.y

  {   Res.z.x :=    ( + S.y.x*S.z.y - S.y.y*S.z.x  );
      Res.z.y :=    ( + S.x.y*S.z.x - S.x.x*S.z.y  );   }

    movhlps XMM10, XMM4
    movlhps XMM10, XMM4    // XMM10 = S.z.y | S.z.x

    pxor XMM12, XMM12
    subpd  XMM12, XMM10   // XMM10 = -S.z.y | -S.z.x

    mulpd   XMM12, XMM0    // XMM12 = S.z.y*S.y.x | S.z.x*S.y.y
    // XMM2 x XMM4
    mulpd XMM10, XMM2     // XMM10 = -S.x.x*S.z.y | -S.x.y*S.z.x
    // XMM0 x XMM4

    hsubpd XMM10, XMM12   // XMM10 = Res.z.x | Res.z.y

  {  Res.x.z :=    ( + S.x.y*S.y.z - S.x.z*S.y.y  ); //  XMM7.L
     Res.y.z :=    ( + S.y.x*S.x.z - S.x.x*S.y.z  ); // XMM9.L

    // Reorder:
    Res.x.z :=    ( + S.x.y*S.y.z - S.y.y*S.x.z  ); //  XMM7.L
    Res.y.z :=   -(   S.x.x*S.y.z - S.y.x*S.x.z  ); // XMM9.L       }

    movlhps XMM3, XMM3   // XMM3 = S.y.z | S.y.z
    mulpd XMM3, XMM0     // XMM3 = S.x.x*S.y.z | S.x.y*S.y.z
    movapd XMM12, XMM2   // XMM12 = S.y.x | S.y.y
    movlhps XMM1, XMM1   // XMM1  = S.x.z | S.x.z

    mulpd XMM12, XMM1    // XMM12 = S.y.x*S.x.z | S.y.y*S.x.z

    subpd XMM3, XMM12    // XMM3 = S.x.x*S.y.z-S.y.x*S.x.z | S.x.y*S.y.z - S.y.y*S.x.z
                         //      = -Res.y.z                | Res.x.z
    pxor XMM9, XMM9      // XMM9 = 0 | 0
    subsd XMM9, XMM3     // XMM9 = Res.y.z
    movhlps XMM7, XMM3   // XMM11 = Res.x.z

  {  Res.z.z :=    ( + S.x.x*S.y.y - S.x.y*S.y.x  );     }

    movhlps XMM11, XMM2    // XMM10 = S.z.y | S.y.y
    movlhps XMM11, XMM2    // XMM10 = S.z.y | S.z.x
    mulpd   XMM11, XMM0    // XMM11 = S.x.x*S.y.y | S.x.y*S.y.x
    hsubpd  XMM11, XMM11   // XMM11 = S.x.x*S.y.y - S.x.y*S.y.x


  {   D := S.x.x*Res.x.x +   S.x.y*Res.y.x +   S.x.z*Res.z.x ==
           Res.x.x*S.x.x +   Res.x.y*S.y.x +   Res.x.z*S.z.x         }

    movlhps XMM0, XMM2
    dppd XMM0, XMM6,  49 // imm8 = bit0 | ~bit1 | bit4 | bit5
                         // bit0- save res to low half
                         // bit4- * low parts
                         // bit5- * high parts

    mulsd XMM4, XMM7     // S.z.x * Res.x.z
    addsd XMM4, XMM0


     //  xmm0 = D

    movsd XMM0, XMM4     // if ABS (D) > eps then PROCEED else EXIT
    movsd XMM1, sign
    pand XMM0, XMM1
    movsd XMM1, eps
    comisd XMM1, XMM0

    jnbe @exit


  {  if Abs(d) < 1.0e-10 then
      Exit;

    D := 1.0 / D;

    S := mult(Res,D);
  }

    movsd XMM3, One
    divsd XMM3, XMM4

    movlhps XMM3, XMM3

    mulpd XMM6, XMM3
    movupd  oWORD[RCX              ], XMM6   //   S.x.x | S.x.y
    mulpd XMM7, XMM3
    movupd  oWORD[RCX+  2*FloatSize], XMM7   //   S.x.z |
    mulpd XMM8, XMM3
    movupd  oWORD[RCX+  4*FloatSize], XMM8   //   S.y.x | S.y.y
    mulpd XMM9, XMM3
    movupd  oWORD[RCX+  6*FloatSize], XMM9   //   S.y.z |
    mulpd XMM10, XMM3
    movupd  oWORD[RCX+  8*FloatSize], XMM10  //   S.z.x | S.z.y
    mulpd XMM11, XMM3
    movupd  oWORD[RCX+ 10*FloatSize], XMM11  //   S.z.z |
  //  sfence
    @exit:
{$ELSE}
begin
  Invert2(S);
{$IFEND}

end;


class procedure T_SSE.Invert(const S: T_Tens; const i0, i1: integer);
{$IFDEF CPUX64}

  const FloatSize = 8;
        One: double = 1.0;
        Zero: double = 0.0;
        eps: double = 1.0e-10;
        sign: Uint64 = $7FFFFFFFFFFFFFFF;
  asm
    .NOFRAME
    .SAVENV XMM5   // S.z.z
    .SAVENV XMM6   // Res.x.x | Res.x.y
    .SAVENV XMM7   // Res.x.z | w
    .SAVENV XMM8   // Res.y.x | Res.y.y
    .SAVENV XMM9   // Res.y.z | w
    .SAVENV XMM10  // Res.z.x | Res.z.y
    .SAVENV XMM11  // Res.z.z | w
    .SAVENV XMM12  // tmp1
    .SAVENV XMM13  // tmp2        }

    XOR RAX, RAX
    mov EAX, i0
    imul RAX, 12*FloatSize
    ADD RCX, RAX // RCX = @S[i0]

    XOR RAX, RAX
    mov EAX, i1
    sub EAX, i0
    inc RAX
    nop


    PXOR XMM6, XMM6
    PXOR XMM7, XMM7
    PXOR XMM8, XMM8
    PXOR XMM9, XMM9
    PXOR XMM10, XMM10
    PXOR XMM11, XMM11

@begin:
    prefetch [RCX+ 12*FloatSize  ]
    prefetch [RCX+ 16*FloatSize  ]
    prefetch [RCX+ 20*FloatSize  ]

    movupd XMM0, oWORD[RCX              ]   //   S.x.x | S.x.y
    movupd XMM1, oWORD[RCX+  2*FloatSize]   //   S.x.z | w
    movupd XMM2, oWORD[RCX+  4*FloatSize]   //   S.y.x | S.y.y
    movupd XMM3, oWORD[RCX+  6*FloatSize]   //   S.y.z | w
    movupd XMM4, oWORD[RCX+  8*FloatSize]   //   S.z.x | S.z.y
    movupd XMM5, oWORD[RCX+ 10*FloatSize]   //   S.z.z | w

  {   Res.x.x :=    ( + S.y.y*S.z.z - S.y.z*S.z.y  );   XMM6
      Res.x.y :=    ( + S.x.z*S.z.y - S.x.y*S.z.z  );          }

    movhlps XMM6, XMM2  // XMM6.L = S.y.y
    movlhps XMM6, XMM3  // XMM6.h = S.y.z
    blendpd XMM5, XMM4, 2  // XMM5= S.z.z | S.z.y

    mulpd XMM6, XMM5   // XMM6 = S.y.y*S.z.z |  S.y.z*S.z.y

    movhlps XMM12, XMM4  // XMM12.L = S.z.y
    movlhps XMM12, XMM5  // XMM12.h = S.z.z
    blendpd XMM1,  XMM0, 2 // XMM1 = S.x.z | S.x.y

    mulpd XMM12, XMM1   // XMM13 = S.x.z*S.z.y | S.x.y* S.z.z

    hsubpd XMM6, XMM12   // XMM6 = Res.x.x | Res.x.y     }

  {  Res.y.x :=    ( + S.z.x*S.y.z - S.y.x*S.z.z  );
     Res.y.y :=    ( + S.x.x*S.z.z - S.x.z*S.z.x  );   }

    movsd   XMM8, XMM4  // XMM8.L = S.z.x
    movlhps XMM8, XMM2  // XMM8.h = S.y.x

    movlhps XMM3 , XMM5  // XMM3 = S.y.z | S.z.z

    mulpd   XMM8, XMM3   // XMM8 = S.z.x*S.y.z |  S.y.x*S.z.z

    movsd   XMM13, XMM0  // XMM12.L = S.x.x
    movlhps XMM13, XMM1  // XMM13 = S.x.x | S.x.z

    movlhps XMM5, XMM4   // XMM5 = S.z.z  | S.z.x

    mulpd XMM13, XMM5    // XMM13 = S.x.x*S.z.z | S.x.z* S.z.x

    hsubpd XMM8, XMM13   // XMM8 = Res.y.x | Res.y.y

  {   Res.z.x :=    ( + S.y.x*S.z.y - S.y.y*S.z.x  );
      Res.z.y :=    ( + S.x.y*S.z.x - S.x.x*S.z.y  );   }

    movhlps XMM10, XMM4
    movlhps XMM10, XMM4    // XMM10 = S.z.y | S.z.x

    pxor XMM12, XMM12
    subpd  XMM12, XMM10   // XMM10 = -S.z.y | -S.z.x

    mulpd   XMM12, XMM0    // XMM12 = S.z.y*S.y.x | S.z.x*S.y.y
    mulpd XMM10, XMM2     // XMM10 = -S.x.x*S.z.y | -S.x.y*S.z.x

    hsubpd XMM10, XMM12   // XMM10 = Res.z.x | Res.z.y

  {  Res.x.z :=    ( + S.x.y*S.y.z - S.x.z*S.y.y  ); //  XMM7.L
     Res.y.z :=    ( + S.y.x*S.x.z - S.x.x*S.y.z  ); // XMM9.L

    // Reorder:
    Res.x.z :=    ( + S.x.y*S.y.z - S.y.y*S.x.z  ); //  XMM7.L
    Res.y.z :=   -(   S.x.x*S.y.z - S.y.x*S.x.z  ); // XMM9.L       }

    movlhps XMM3, XMM3   // XMM3 = S.y.z | S.y.z
    mulpd XMM3, XMM0     // XMM3 = S.x.x*S.y.z | S.x.y*S.y.z

    movapd XMM12, XMM2   // XMM12 = S.y.x | S.y.y
    movlhps XMM1, XMM1   // XMM1  = S.x.z | S.x.z

    mulpd XMM12, XMM1    // XMM12 = S.y.x*S.x.z | S.y.y*S.x.z

    subpd XMM3, XMM12    // XMM3 = S.x.x*S.y.z-S.y.x*S.x.z | S.x.y*S.y.z - S.y.y*S.x.z
                         //      = -Res.y.z                | Res.x.z
    pxor XMM9, XMM9      // XMM9 = 0 | 0
    subsd XMM9, XMM3     // XMM9 = Res.y.z
    movhlps XMM7, XMM3   // XMM11 = Res.x.z

  {  Res.z.z :=    ( + S.x.x*S.y.y - S.x.y*S.y.x  );     }

    movhlps XMM11, XMM2    // XMM10 = S.y.y | w
    movlhps XMM11, XMM2    // XMM10 = S.y.y | S.y.x
    mulpd   XMM11, XMM0    // XMM11 = S.x.x*S.y.y | S.x.y*S.y.x
    hsubpd  XMM11, XMM11   // XMM11 = S.x.x*S.y.y - S.x.y*S.y.x


  {   D := S.x.x*Res.x.x +   S.x.y*Res.y.x +   S.x.z*Res.z.x ==
           Res.x.x*S.x.x +   Res.x.y*S.y.x +   Res.x.z*S.z.x         }

    movlhps XMM0, XMM2
    dppd XMM0, XMM6,  49 // imm8 = bit0 | ~bit1 | bit4 | bit5
                         // bit0- save res to low half
                         // bit4- * low parts
                         // bit5- * high parts

    mulsd XMM4, XMM7     // S.z.x * Res.x.z
    addsd XMM4, XMM0

     //  xmm4 = D

    movsd XMM0, XMM4     // if ABS (D) > eps then PROCEED else EXIT
    movsd XMM1, sign
    pand XMM0, XMM1
    movsd XMM1, eps
    comisd XMM1, XMM0

    jnbe @exit

  //  D := 1.0 / D;
  //  S := mult(Res,D);

    movsd XMM3, One
    divsd XMM3, XMM4

    movlhps XMM3, XMM3

    mulpd XMM6, XMM3
    movupd  oWORD[RCX              ], XMM6   //   S.x.x | S.x.y
    add RCX, 2*FloatSize

    mulsd XMM7, XMM3
    movupd  [RCX        ], XMM7   //   S.x.z |
    add RCX, 2*FloatSize

    mulpd XMM8, XMM3
    movupd  [RCX         ], XMM8   //   S.y.x | S.y.y
    add RCX, 2*FloatSize

    mulsd XMM9, XMM3
    movupd   [RCX    ], XMM9   //   S.y.z |
    add RCX, 2*FloatSize

    mulpd XMM10, XMM3
    movupd  [RCX      ], XMM10  //   S.z.x | S.z.y
    add RCX, 2*FloatSize

    mulsd XMM11, XMM3
    movupd   [RCX     ], XMM11  //   S.z.z |
    add RCX, 2*FloatSize
  //  sfence
    @exit:

//    Add RCX, 12*FloatSize  // next matrix
    dec RAX

    cmp RAX, 0
    jg @begin


{$ELSE}
var i: integer;
begin
  for i := i0 to i1 do
    Invert2(  T_TensorArr(@S)[i] );
{$IFEND}
end;


class procedure T_SSE.Invert_gauss(var S: T_M4);
{$IFDEF CPUX64}
const FloatSize = 8;
      One: double = 1.0;
      Zero: double = 0.0;
asm

  .NOFRAME
  .SAVENV XMM5
  .SAVENV XMM6
  .SAVENV XMM7
  .SAVENV XMM8
  .SAVENV XMM9
  .SAVENV XMM10
//  .SAVENV XMM11   // }


  movupd XMM0, [RCX              ]   //   S.x.x | S.x.y
  movupd XMM1, [RCX+  2*FloatSize]   //   S.x.z | S.x.t

  prefetch [RCX + 4*FloatSize]     // slightly increase overall performance (~1%)

  PXOR XMM10,XMM10

//****************  FORWARD MOVE ***********************

 // STEP  #1
 // STEP #1.prep:  D = 1.0/ S.x.x    S.x.x = 1; S.x := S.x*D;

  movsd XMM8, One             // XMM8.L = 1.0;
  movsd XMM9, XMM0          // EXTRACT DIVIDER!
  movsd XMM0, XMM8            //XMM0 = 1 | S.x.y
  divsd XMM8, XMM9            // D := 1.0 / S.x.x;
  movlhps XMM8, XMM8          // XMM8 = D  | D


  mulpd XMM0, XMM8            // S.x = S.x*D
  mulpd XMM1, XMM8            // XMM0 XMM1= a00*d | a01*d   ||  a02*d | a03*d

  movupd XMM2, [RCX+  4*FloatSize]   //   S.y.x | S.y.y
  movupd XMM3, [RCX+  6*FloatSize]   //   S.y.z | S.y.t

  prefetch [RCX + 8*FloatSize]



 // STEP  #1.y:   D := S.y.x;  s.y.x :=  0; S.y = S.y- S.x*D

  movsd   XMM8, XMM2   // XMM8 = (D= S.y.x) | 0
  movlhps XMM8, XMM8   // XMM8 =  D          |   D

  movsd XMM2, XMM10    // S.y.x := 0;  XMM2 = 0 | S.y.y

  movapd XMM9, XMM0    // copy S.x.x, S.x.y
  mulpd XMM9, XMM8     // mul
  subpd XMM2, XMM9     // XMM2 = -S.x.x*D         | S.y.y - S.x.y*D

  movapd XMM9, XMM1    // copy S.x.z, S.x.t
  mulpd XMM9, XMM8     // mul
  subpd XMM3, XMM9     // XMM3 = S.y.z - S.x.z*d |  S.y.t - S.x.t*d        //}


 // STEP  #1.z:   D := S.z.x;  s.z.x :=  0; S.z = S.z- S.x*D
  movupd XMM4, [RCX+  8*FloatSize]   //   S.z.x | S.z.y
  movupd XMM5, [RCX+ 10*FloatSize]   //   S.z.z | S.z.t

  prefetch [RCX + 12*FloatSize]


  movsd   XMM8, XMM4   // XMM8 = (D= S.y.x) | 0
  movlhps XMM8, XMM8   // XMM8 =  D          |   D

  movsd XMM4, XMM10    // S.z.x:= 0;  XMM4 = 0 | S.z.y

  movapd XMM9, XMM0    // Copy S.x.x, S.x.y
  mulpd XMM9, XMM8     // mul
  subpd XMM4, XMM9     // XMM4 = -S.x.x*D         | S.y.y - S.x.y*D

  movapd XMM9, XMM1    // copy S.x.z , S.x.t
  mulpd XMM9, XMM8     // mul
  subpd XMM5, XMM9     // XMM5 = S.z.z - S.z.z*d |  S.z.t - S.z.t*d        //}

 // STEP  #1.t:   D := S.t.x;  s.t.x :=  0; S.t = S.t- S.x*D
  movupd XMM6, [RCX+ 12*FloatSize]   //   S.t.x | S.t.y
  movupd XMM7, [RCX+ 14*FloatSize]   //   S.t.z | S.t.t

  movsd   XMM8, XMM6   // XMM8 = (D= S.y.x) | 0
  movlhps XMM8, XMM8   // XMM8 =  D          |   D

  movsd XMM6, XMM10    // S.t.x:= 0;  XMM6 = 0 | S.z.y

  movapd XMM9, XMM0    // Copy S.y.x, S.y.y
  mulpd XMM9, XMM8     // mul
  subpd XMM6, XMM9     // XMM6 = -S.x.x*D         | S.y.y - S.x.y*D

  movapd XMM9, XMM1    // copy S.z.z , S.z.t
  mulpd XMM9, XMM8     // mul
  subpd XMM7, XMM9     // XMM7 = S.z.z - S.z.z*d |  S.z.t - S.z.t*d        //}


// STEP  #4.prep  D := 1.0/S.y.y;    S.y.y = 1;  S.y := S.y*D;

  movsd   XMM8, One
  movhlps XMM9, XMM2    // XMM9 = S.y.y    |    w
  movlhps XMM2, XMM8    // XMM2 = S.y.x    |    S.y.y=1
  divsd   XMM8, XMM9    // XMM8 = 1.0/ S.y.y   = D
  movlhps XMM8, XMM8

  mulpd   XMM2, XMM8    // XMM2 = S.y.x*D  |    D= 1/S.y.y
  mulpd   XMM3, XMM8    // XMM3 = S.y.z*D  |    S.y.t*D

// STEP  #4.z  D := S.z.y; S.z.y := 0; S.z := S.z-S.y*D

  movhlps XMM8, XMM4   // XMM9= D=S.z.y   |   w

  movlhps XMM8, XMM8   // XMM8= D         |   D

  movlhps XMM4, XMM10  // XMM4= S.z.x     | S.z.y=0

  movapd  XMM9, XMM2   // XMM9= S.y.x      | S.y.y
  mulpd   XMM9, XMM8   // XMM9= S.y.x*D    | S.y.y*D
  subpd   XMM4, XMM9   // XMM4= S.z.x - S.y.x*D | 0 -S.y.y*D

  movapd  XMM9, XMM3   // XMM9= S.y.z      | S.y.t
  mulpd   XMM9, XMM8   // XMM9= S.y.z*D    | S.y.t*D
  subpd   XMM5, XMM9   // XMM4= S.z.z - S.y.z*D | S.z.t - S.y.t*D


// STEP  #4.t  D := S.t.y; S.t.y := 0; S.t := S.t-S.y*D

  movhlps XMM8, XMM6   // XMM9= D=S.t.y   |   w

  movlhps XMM8, XMM8   // XMM8= D         |   D

  movlhps XMM6, XMM10  // XMM4= S.t.x     | S.t.y=0

  movapd  XMM9, XMM2   // XMM9= S.y.x      | S.y.y
  mulpd   XMM9, XMM8   // XMM9= S.y.x*D    | S.y.y*D
  subpd   XMM6, XMM9   // XMM4= S.t.x - S.y.x*D | 0 -S.y.y*D

  movapd  XMM9, XMM3   // XMM9= S.y.z      | S.y.t
  mulpd   XMM9, XMM8   // XMM9= S.y.z*D    | S.y.t*D
  subpd   XMM7, XMM9   // XMM4= S.t.z - S.y.z*D | S.t.t - S.y.t*D

// STEP  #5.prep  D := 1.0/S.z.z;    S.z.z = 1;  S.z := S.z*D;

  movsd   XMM8, One
  movsd   XMM9, XMM5    // XMM9 = S.z.z    |    w
  movsd   XMM5, XMM8    // XMM5 = S.z.z=1  |    S.z.t
  divsd   XMM8, XMM9    // XMM8 = 1.0/ S.z.z   = D
  movlhps XMM8, XMM8

  mulpd   XMM4, XMM8    // XMM4 = S.z.x*D  |    D= 1/S.z.z
  mulpd   XMM5, XMM8    // XMM5 = S.z.z*D  |    S.z.t*D

// STEP  #5.t  D := S.t.z; S.t.z := 0; S.t := S.t-S.z*D

  movsd   XMM8, XMM7   // XMM8= D=S.t.z   |   w

  movlhps XMM8, XMM8   // XMM8= D         |   D

  movsd   XMM7, XMM10  // XMM7= S.t.z=0   | S.t.t

  movapd  XMM9, XMM4   // XMM9= S.z.x      | S.z.y
  mulpd   XMM9, XMM8   // XMM9= S.z.x*D    | S.z.y*D
  subpd   XMM6, XMM9   // XMM4= S.t.x - S.z.x*D | S.t.y -S.z.y*D

  movapd  XMM9, XMM5   // XMM9= S.z.z      | S.z.t
  mulpd   XMM9, XMM8   // XMM9= S.z.z*D    | S.z.t*D
  subpd   XMM7, XMM9   // XMM7= 0 - S.z.z*D | S.t.t - S.z.t*D

// STEP  #6  D := 1/S.t.t; S.t.t := 1.0; S.t := S.t*D

  movsd   XMM8, One
  movhlps XMM9, XMM7    // XMM9 = S.t.t    |    w
  movlhps XMM7, XMM8    // XMM7 = S.t.z    |    S.t.t=1
  divsd   XMM8, XMM9    // XMM8 = 1.0/ S.z.z   = D
  movlhps XMM8, XMM8

  mulpd XMM6, XMM8
  mulpd XMM7, XMM8


  // save S.t
  movupd oWORD[RCX+ 12*FloatSize], XMM6   //   S.z.x | S.z.y
  movupd oWORD[RCX+ 14*FloatSize], XMM7   //   S.z.z |


// BACKWARD MOVEMENT!

// STEP  #7.x-y  D := S.x.y; S.x.y := 0; S.x -= S.y*D

  movhlps XMM8, XMM0     // XMM8= D=S.x.y | w
  movlhps XMM8, XMM8     // XMM8 = D | D

  movlhps XMM0, XMM10    // S.x.y = 0
  movapd XMM9, XMM2      // XMM9= S.y.x      | S.y.y
  mulpd XMM9, XMM8       // XMM9= S.y.x*D    | S.y.y*D
  subpd XMM0, XMM9       // XMM0= S.x.x+S.y.x*D | S.y.y*D

  movapd XMM9, XMM3       // XMM9= S.y.z      | t
  mulpd XMM9, XMM8       // XMM9= S.y.z*D    | t
  subpd XMM1, XMM9       // XMM0= S.x.z+S.y.z*D | t

// STEP  #7.y-z  D := S.y.z; S.y.z := 0; S.y -= S.z*D

  movsd   XMM8, XMM3     // XMM8= D=S.y.z | w
  movsd   XMM3, XMM10    // S.y.z := 0
  movlhps XMM8, XMM8     // XMM8 = D | D


  movapd XMM9, XMM4      // XMM9= S.z.x      | S.z.y
  mulpd XMM9, XMM8       // XMM9= S.z.x*D    | S.z.y*D
  subpd XMM2, XMM9       // XMM0= S.y.x-S.z.x*D | S.y.y-S.z.y*D

  movapd XMM9, XMM5       // XMM9= S.y.z      | t
  mulpd XMM9, XMM8       // XMM9= S.y.z*D    | t
  subpd XMM3, XMM9       // XMM0= S.x.z+S.y.z*D | t

// STEP  #7.x-z  D := S.x.z; S.x.z := 0; S.x -= S.z*D

  movsd   XMM8, XMM1     // XMM8= D=S.x.z | w
  movsd   XMM1, XMM10    // S.x.z := 0
  movlhps XMM8, XMM8     // XMM8 = D | D


  movapd XMM9, XMM4      // XMM9= S.t.x      | S.t.y
  mulpd XMM9, XMM8       // XMM9= S.t.x*D    | S.t.y*D
  subpd XMM0, XMM9       // XMM0= S.z.x-S.t.x*D | S.z.y-S.t.y*D

  movapd XMM9, XMM5      // XMM9= S.t.z      | S.t.t
  mulpd XMM9, XMM8       // XMM9= S.t.z*D    | S.t.t*D
  subpd XMM1, XMM9       // XMM0= S.z.z-S.t.z*D | S.z.t-S.t.t*D

// STEP  #7.z-t  D := S.z.t; S.z.t := 0; S.z -= S.t*D

  movhlps XMM8, XMM5     // XMM8= D=S.z.t | w
  movlhps XMM5, XMM10    // S.z.t := 0
  movlhps XMM8, XMM8     // XMM8 = D | D


  movapd XMM9, XMM6      // XMM9= S.t.x      | S.t.y
  mulpd XMM9, XMM8       // XMM9= S.t.x*D    | S.t.y*D
  subpd XMM4, XMM9       // XMM0= S.z.x-S.t.x*D | S.z.y-S.t.y*D

  movapd XMM9, XMM7      // XMM9= S.t.z      | S.t.t
  mulpd XMM9, XMM8       // XMM9= S.t.z*D    | S.t.t*D
  subpd XMM5, XMM9       // XMM0= S.z.z-S.t.z*D | S.z.t-S.t.t*D

  // SAVE S.z
  movupd oWORD[RCX+  8*FloatSize], XMM4   //   S.z.x | S.z.y
  movupd oWORD[RCX+ 10*FloatSize], XMM5   //   S.z.z |
// STEP  #8.y-t  D := S.y.t; S.y.t := 0; S.y -= S.t*D

  movhlps XMM8, XMM3     // XMM8= D=S.y.t | w
  movlhps XMM3, XMM10    // S.y.t = 0

  movlhps XMM8, XMM8     // XMM8 = D | D

  movapd XMM9, XMM6      // XMM9= S.y.x      | S.y.y
  mulpd XMM9, XMM8       // XMM9= S.y.x*D    | S.y.y*D
  subpd XMM2, XMM9       // XMM0= S.x.x+S.y.x*D | S.y.y*D

  movapd XMM9, XMM7      // XMM9= S.y.z      | t
  mulpd XMM9, XMM8       // XMM9= S.y.z*D    | t
  subpd XMM3, XMM9       // XMM0= S.x.z+S.y.z*D | t

//   save S.y
  movupd oWORD[RCX+  4*FloatSize], XMM2   //   S.y.x | S.y.y
  movupd oWORD[RCX+  6*FloatSize], XMM3   //   S.y.z |

// STEP  #8.x-t  D := S.x.t; S.x.t := 0; S.x -= S.t*D

  movhlps XMM8, XMM1     // XMM8= D=S.x.t | w
  movlhps XMM1, XMM10    // S.x.t := 0
  movlhps XMM8, XMM8     // XMM8 = D | D


  movapd XMM9, XMM6      // XMM9= S.t.x      | S.t.y
  mulpd XMM9, XMM8       // XMM9= S.t.x*D    | S.t.y*D
  subpd XMM0, XMM9       // XMM0= S.x.x-S.t.x*D | S.x.y-S.t.y*D

  movapd XMM9, XMM7      // XMM9= S.t.z         | S.t.t
  mulpd XMM9, XMM8       // XMM9= S.t.z*D       | S.t.t*D
  subpd XMM1, XMM9       // XMM0= S.x.z-S.t.z*D | S.x.t-S.t.t*D



  movupd oWORD[RCX              ], XMM0   //   S.x.x | S.x.y
  movupd oWORD[RCX+  2*FloatSize], XMM1   //   S.x.z |
 {$ELSE}
begin
  Invert_m4( S );
 {$ENDIF}
end;


class procedure T_SSE.Invert_gauss(var S: T_Tens);
{$ifdef CPUX64}
const FloatSize = 8;
      One: double = 1.0;
      Zero: double = 0.0;
asm
// STRAIGTH MOVEMENT!
  // STEP  #1
  .NOFRAME
  .SAVENV XMM5
  .SAVENV XMM6
  .SAVENV XMM7
  .SAVENV XMM8
  .SAVENV XMM9
  .SAVENV XMM10
  .SAVENV XMM11   // }


  movapd XMM0, oWORD ptr [RCX              ]   //   S.x.x | S.x.y
  movapd XMM1, oWORD ptr [RCX+  2*FloatSize]   //   S.x.z | w
  movapd XMM2, oWORD ptr [RCX+  4*FloatSize]   //   S.y.x | S.y.y
  movapd XMM3, oWORD ptr [RCX+  6*FloatSize]   //   S.y.z | w
  movapd XMM4, oWORD ptr [RCX+  8*FloatSize]   //   S.z.x | S.z.y
  movapd XMM5, oWORD ptr [RCX+ 10*FloatSize]   //   S.z.z | w

//  PXOR XMM6, XMM6

{  // STEP #1

  D = 1.0/ S.x.x
  S.x.x :=       D;
  S.x.y := S.x.y*D;
  S.x.z := S.x.z*D;}



  movsd XMM11, One
      //  movlhps XMM11, XMM11   //  XMM11 = 1 | 1
      //  movsd XMM10, Zero
      //  movlhps XMM10, XMM10   //  XMM10 = 0 | 0
  PXOR XMM10,XMM10

  movsd XMM6, XMM11                 // XMM6.L = 1.0;
  divsd XMM6, XMM0                  // D := 1.0 / S.x.x;
  movlhps XMM6, XMM6                // XMM6 = D  | D

  mulpd XMM0, XMM6
  mulsd XMM1, XMM6                  // XMM0 XMM1= a00*d | a01*d   ||  a02*d | x

  movsd XMM0, XMM6                  // XMM0, XMM1 = D, a01*d, a02*d, x


{    // STEP  #2
  D := S.y.x;
  s.y.x :=  0    - S.x.x*d;
  s.y.y := S.y.y - S.x.y*d;
  s.y.z := S.y.z - S.x.z*d;}

{  PXOR XMM6,XMM6     // XMM6 = 0 | 0
  subsd   XMM6, XMM2   // XMM6 = (D= -S.y.x) | 0
  movlhps XMM6, XMM6   // XMM6 =  D          |   D

  movsd XMM2, XMM10     // XMM2 = 0 | S.y.y

  movapd XMM7, XMM0
  mulpd XMM7, XMM6

  addpd XMM2, XMM7     // XMM2 = S.x.x*D         | S.y.y + S.x.y*D

  movsd XMM7, XMM1
  mulsd XMM7, XMM6
  addsd XMM3, XMM7     // XMM3 = S.y.z + S.x.z*d |  w    // }



  movsd   XMM6, XMM2   // XMM6 = (D= S.y.x) | 0
  movlhps XMM6, XMM2   // XMM6 =  D          |   D

  movsd XMM2, XMM10     // XMM2 = 0 | S.y.y

  movapd XMM7, XMM0
  mulpd XMM7, XMM6

  subpd XMM2, XMM7     // XMM2 = -S.x.x*D         | S.y.y - S.x.y*D

  movsd XMM7, XMM1
  mulsd XMM7, XMM6
  subsd XMM3, XMM7     // XMM3 = S.y.z - S.x.z*d |  w        //}

{    // STEP  #3
  D := -S.z.x;
  S.z.x :=         S.x.x*d;
  S.z.y := S.z.y + S.x.y*d;
  S.z.z := S.z.z + S.x.z*d;}

//  PXOR XMM6,XMM6       // XMM6 = 0 | 0

  movsd   XMM6, XMM4   // XMM6 = (D=  S.z.x) | 0
  movlhps XMM6, XMM4   // XMM6 =  D          |   D

  movsd XMM4, XMM10    // XMM4 = 0 | S.z.y

  movapd XMM7, XMM0
  mulpd XMM7, XMM6

  subpd XMM4, XMM7     // XMM4 = -S.z.x*D         | S.z.y - S.x.y*D

  movsd XMM7, XMM1
  mulsd XMM7, XMM6
  subsd XMM5, XMM7     // XMM5 = S.z.z - S.x.z*d |  w


{    // STEP  #4
  D := 1.0/S.y.y;
  S.y.x := S.y.x*D;
  S.y.y :=       D;
  S.y.z := S.y.z*D;  }

  movsd   XMM6, XMM11
  movhlps XMM7, XMM2    // XMM7 = S.y.y    |    w
  divsd   XMM6, XMM7    // XMM6 = 1.0/ S.y.y   = D
  movlhps XMM6, XMM6

  mulsd   XMM2, XMM6    // XMM2 = S.y.x*D  |   S.y.y
  movlhps XMM2, XMM6    // XMM2 = S.y.x*D  |   D

  mulsd   XMM3, XMM6    // XMM3 = S.y.z*D  |    w

{    // STEP  #5
  D := -S.z.y;                  //  S.z = XMM4 | XMM5
  S.z.x := S.z.x + S.y.x*D;
  S.z.y :=         S.y.y*D;
  S.z.z := S.z.z + S.y.z*D;}


  movhlps XMM6, XMM4   // XMM7= S.z.y     |   w
  movlhps XMM6, XMM6   // XMM6= D         |   D

  movlhps XMM4, XMM10  // XMM4= S.z.x     | 0

  movapd   XMM7, XMM2   // XMM7= S.y.x      | S.y.y
  mulpd   XMM7, XMM6   // XMM7= S.y.x*D    | S.y.y*D
  subpd   XMM4, XMM7   // XMM4= S.z.x - S.y.x*D | -S.y.y*D

//  movsd   XMM7, XMM3   // XMM7= S.y.z      | t
//  mulsd   XMM7, XMM6   // XMM7= S.y.z*D    | t
//  subsd   XMM5, XMM7   // XMM5= S.z.z - S.y.z*D | t

//  movsd   XMM7, XMM3   // XMM7= S.y.z      | t
  mulsd   XMM6, XMM3   // XMM6= S.y.z*D    | t
  subsd   XMM5, XMM6   // XMM5= S.z.z - S.y.z*D | t

 {   // STEP  #6
  D := 1.0/S.z.z;
  S.z.x := S.z.x*D;
  S.z.y := S.z.y*D;
  S.z.z :=       D;    // 1/A22}

//  movsd XMM7, XMM11    // XMM7= 1 | w
  divsd XMM11, XMM5     // XMM11= 1/S.z.z  | 1.0
  movsd XMM5, XMM11
  movlhps XMM5, XMM5   // XMM7 = D  | D

  mulpd XMM4, XMM5     // XMM4 = S.z.x*D | S.z.y*D

// BACKWARD MOVEMENT!

{    // STEP  #7
  D := -S.x.y;
  S.x.x := S.x.x + S.y.x*D;
  S.x.y :=         S.y.y*D;
  S.x.z := S.x.z + S.y.z*D;}

  movapd XMM6, XMM0      // XMM7= S.x.x | S.x.y
  movhlps XMM6, XMM6     // XMM6 = D | D   =S.x.y

  movlhps XMM0, XMM10    // S.x.y = 0
  movapd XMM7, XMM2      // XMM7= S.y.x      | S.y.y
  mulpd XMM7, XMM6       // XMM7= S.y.x*D    | S.y.y*D
  subpd XMM0, XMM7       // XMM0= S.x.x+S.y.x*D | S.y.y*D

  movsd XMM7, XMM3       // XMM7= S.y.z      | t
  mulsd XMM7, XMM6       // XMM7= S.y.z*D    | t
  subsd XMM1, XMM7       // XMM0= S.x.z+S.y.z*D | t


{    // STEP  #8
  D := S.y.z;
  S.y.x := S.y.x - S.z.x*D;     Xmm2.L
  S.y.y := S.y.y - S.z.y*D;     Xmm2.H
  S.y.z :=       - S.z.z*D;     Xmm3.L
  }

//  movsd XMM6, XMM3       // XMM7= S.y.z= D
  movlhps XMM3, XMM3     // XMM3 = D | D

//  movlhps XMM0, XMM10  // S.x.y = 0
  movapd XMM7, XMM4      // XMM7= S.z.x         | S.z.y
  mulpd XMM7, XMM3       // XMM7= S.z.x*D       | S.z.y*D
  subpd XMM2, XMM7       // XMM2= S.y.x-S.z.x*D | S.y.y -  S.z.y*D

  movsd XMM7, XMM5       // XMM7= S.z.z          | t
  mulsd XMM7, XMM3       // XMM7= S.z.z*D
  movsd XMM3, XMM10
  subsd XMM3, XMM7


{    // STEP  #9
  D := -S.x.z;
  S.x.x := S.x.x + S.z.x*D;
  S.x.y := S.x.y + S.z.y*D;
  S.x.z :=         S.z.z*D;}

  movsd XMM6, XMM10    // 0
  subsd XMM6, XMM1     // 0  -S.x.z == D
  movlhps XMM6, XMM6   // XMM6 = D  |  D

  movapd XMM7, XMM4    // S.z.x    | S.z.y
  mulpd  XMM7, XMM6    // S.z.x*D  | S.z.y*D

  addpd XMM0, XMM7     // XMM0= S.x.x + S.z.x*D  | S.x.y + S.z.y*D
  movsd XMM1, XMM5     // XMM1= S.z.z
  mulsd XMM1, XMM6     // XMM1= S.z.z*D


  movapd  oWORD ptr[RCX              ], XMM0   //   S.x.x | S.x.y
  movapd  oWORD ptr[RCX+  2*FloatSize], XMM1   //   S.x.z |
  movapd  oWORD ptr[RCX+  4*FloatSize], XMM2   //   S.y.x | S.y.y
  movapd  oWORD ptr[RCX+  6*FloatSize], XMM3   //   S.y.z |
  movapd  oWORD ptr[RCX+  8*FloatSize], XMM4   //   S.z.x | S.z.y
  movapd  oWORD ptr[RCX+ 10*FloatSize], XMM5   //   S.z.z |

{$ELSE}
begin
  Invert2(  S );
{$ENDIF CPUX64}

end;


class procedure T_SSE.Invert_gauss(var S: T_M4; const N: integer);
{$IFDEF CPUX64}
const FloatSize = 8;
      One: double = 1.0;
      Zero: double = 0.0;
asm

  .NOFRAME
  .SAVENV XMM5
  .SAVENV XMM6
  .SAVENV XMM7
  .SAVENV XMM8
  .SAVENV XMM9
  .SAVENV XMM10
  .PUSHNV RBX

  XOR RBX, RBX
  mov EBX, N

@next:


//************** FORWARD MOVE ****************

  // STEP  #1

  MOVNTDQA  XMM0, [RCX              ]   //   S.x.x | S.x.y
  MOVNTDQA  XMM1, [RCX+  2*FloatSize]   //   S.x.z | S.x.t

  prefetch [RCX + 4*FloatSize]
  PXOR XMM10,XMM10

 // STEP #1.prep:  D = 1.0/ S.x.x    S.x.x = 1; S.x := S.x*D;

           //  movsd XMM8, One             // XMM8.L = 1.0;
           //  movsd XMM0, XMM8          //XMM0 = 1 | S.x.y
           //  divsd XMM8, XMM9          // D := 1.0 / S.x.x;
           //  movlhps XMM8, XMM8          // XMM8 = D  | D
           //  mulpd XMM0, XMM8            // S.x = S.x*D
           //  mulpd XMM1, XMM8            // XMM0 XMM1= a00*d | a01*d   ||  a02*d | a03*d

  movapd XMM9, XMM0          // EXTRACT DIVIDER!
  movlhps XMM9, XMM9
  movlpd XMM0, One           // XMM0 = 1.0 | a01


  divpd XMM0, XMM9          //  XMM0 = 1.0/a00 | a01/a00

  movapd XMM8, XMM0         //  XMM8 = 1/a00 =d
  movlhps XMM8, XMM8        // xmm8 = d | d

  mulpd XMM1, XMM8            // XMM0 XMM1= 1.0/a00 | a01*d   ||  a02*d | a03*d

  MOVNTDQA  XMM2, [RCX+  4*FloatSize]   //   S.y.x | S.y.y
  MOVNTDQA  XMM3, [RCX+  6*FloatSize]   //   S.y.z | S.y.t

  prefetch [RCX + 8*FloatSize]

 // STEP  #1.y:   D := S.y.x;  s.y.x :=  0; S.y = S.y- S.x*D

  movapd   XMM8, XMM2   // XMM8 = (D= S.y.x) | 0
  movlhps XMM8, XMM8   // XMM8 =  D          |   D

  movhlps XMM2, XMM10    // S.y.x := 0;  XMM2 = 0 | S.y.y

  movapd XMM9, XMM0    // copy S.x.x, S.x.y
  mulpd XMM9, XMM8     // mul
  subpd XMM2, XMM9     // XMM2 = -S.x.x*D         | S.y.y - S.x.y*D

  movapd XMM9, XMM1    // copy S.x.z, S.x.t
  mulpd XMM9, XMM8     // mul
  subpd XMM3, XMM9     // XMM3 = S.y.z - S.x.z*d |  S.y.t - S.x.t*d        //}


 // STEP  #1.z:   D := S.z.x;  s.z.x :=  0; S.z = S.z- S.x*D
  MOVNTDQA  XMM4, [RCX+  8*FloatSize]   //   S.z.x | S.z.y
  MOVNTDQA  XMM5, [RCX+ 10*FloatSize]   //   S.z.z | S.z.t

  prefetch [RCX + 12*FloatSize]


  movapd   XMM8, XMM4   // XMM8 = (D= S.y.x) | 0
  movlhps XMM8, XMM8   // XMM8 =  D          |   D

  movsd XMM4, XMM10    // S.z.x:= 0;  XMM4 = 0 | S.z.y

  movapd XMM9, XMM0    // Copy S.x.x, S.x.y
  mulpd XMM9, XMM8     // mul
  subpd XMM4, XMM9     // XMM4 = -S.x.x*D         | S.y.y - S.x.y*D

  movapd XMM9, XMM1    // copy S.x.z , S.x.t
  mulpd XMM9, XMM8     // mul
  subpd XMM5, XMM9     // XMM5 = S.z.z - S.z.z*d |  S.z.t - S.z.t*d        //}

 // STEP  #1.t:   D := S.t.x;  s.t.x :=  0; S.t = S.t- S.x*D
  MOVNTDQA  XMM6, [RCX+ 12*FloatSize]   //   S.t.x | S.t.y
  MOVNTDQA  XMM7, [RCX+ 14*FloatSize]   //   S.t.z | S.t.t

  movapd   XMM8, XMM6   // XMM8 = (D= S.y.x) | 0
  movlhps XMM8, XMM8   // XMM8 =  D          |   D

  movhlps XMM6, XMM10    // S.t.x:= 0;  XMM6 = 0 | S.z.y

  movapd XMM9, XMM0    // Copy S.y.x, S.y.y
  mulpd XMM9, XMM8     // mul
  subpd XMM6, XMM9     // XMM6 = -S.x.x*D         | S.y.y - S.x.y*D

  movapd XMM9, XMM1    // copy S.z.z , S.z.t
  mulpd XMM9, XMM8     // mul
  subpd XMM7, XMM9     // XMM7 = S.z.z - S.z.z*d |  S.z.t - S.z.t*d        //}


// STEP  #4.prep  D := 1.0/S.y.y;    S.y.y = 1;  S.y := S.y*D;

  movsd   XMM8, One
  movhlps XMM9, XMM2    // XMM9 = S.y.y    |    w
  movlhps XMM2, XMM8    // XMM2 = S.y.x    |    S.y.y=1
  divsd   XMM8, XMM9    // XMM8 = 1.0/ S.y.y   = D
  movlhps XMM8, XMM8

  mulpd   XMM2, XMM8    // XMM2 = S.y.x*D  |    D= 1/S.y.y
  mulpd   XMM3, XMM8    // XMM3 = S.y.z*D  |    S.y.t*D

// STEP  #4.z  D := S.z.y; S.z.y := 0; S.z := S.z-S.y*D

  movhlps XMM8, XMM4   // XMM9= D=S.z.y   |   w

  movlhps XMM8, XMM8   // XMM8= D         |   D

  movlhps XMM4, XMM10  // XMM4= S.z.x     | S.z.y=0

  movapd  XMM9, XMM2   // XMM9= S.y.x      | S.y.y
  mulpd   XMM9, XMM8   // XMM9= S.y.x*D    | S.y.y*D
  subpd   XMM4, XMM9   // XMM4= S.z.x - S.y.x*D | 0 -S.y.y*D

  movapd  XMM9, XMM3   // XMM9= S.y.z      | S.y.t
  mulpd   XMM9, XMM8   // XMM9= S.y.z*D    | S.y.t*D
  subpd   XMM5, XMM9   // XMM4= S.z.z - S.y.z*D | S.z.t - S.y.t*D


// STEP  #4.t  D := S.t.y; S.t.y := 0; S.t := S.t-S.y*D

  movhlps XMM8, XMM6   // XMM9= D=S.t.y   |   w

  movlhps XMM8, XMM8   // XMM8= D         |   D

  movlhps XMM6, XMM10  // XMM4= S.t.x     | S.t.y=0

  movapd  XMM9, XMM2   // XMM9= S.y.x      | S.y.y
  mulpd   XMM9, XMM8   // XMM9= S.y.x*D    | S.y.y*D
  subpd   XMM6, XMM9   // XMM4= S.t.x - S.y.x*D | 0 -S.y.y*D

  movapd  XMM9, XMM3   // XMM9= S.y.z      | S.y.t
  mulpd   XMM9, XMM8   // XMM9= S.y.z*D    | S.y.t*D
  subpd   XMM7, XMM9   // XMM4= S.t.z - S.y.z*D | S.t.t - S.y.t*D

// STEP  #5.prep  D := 1.0/S.z.z;    S.z.z = 1;  S.z := S.z*D;

  movsd   XMM8, One
  movsd   XMM9, XMM5    // XMM9 = S.z.z    |    w
  movsd   XMM5, XMM8    // XMM5 = S.z.z=1  |    S.z.t
  divsd   XMM8, XMM9    // XMM8 = 1.0/ S.z.z   = D
  movlhps XMM8, XMM8

  mulpd   XMM4, XMM8    // XMM4 = S.z.x*D  |    D= 1/S.z.z
  mulpd   XMM5, XMM8    // XMM5 = S.z.z*D  |    S.z.t*D

// STEP  #5.t  D := S.t.z; S.t.z := 0; S.t := S.t-S.z*D

  movsd   XMM8, XMM7   // XMM8= D=S.t.z   |   w

  movlhps XMM8, XMM8   // XMM8= D         |   D

  movsd   XMM7, XMM10  // XMM7= S.t.z=0   | S.t.t

  movapd  XMM9, XMM4   // XMM9= S.z.x      | S.z.y
  mulpd   XMM9, XMM8   // XMM9= S.z.x*D    | S.z.y*D
  subpd   XMM6, XMM9   // XMM4= S.t.x - S.z.x*D | S.t.y -S.z.y*D

  movapd  XMM9, XMM5   // XMM9= S.z.z      | S.z.t
  mulpd   XMM9, XMM8   // XMM9= S.z.z*D    | S.z.t*D
  subpd   XMM7, XMM9   // XMM7= 0 - S.z.z*D | S.t.t - S.z.t*D

// STEP  #6  D := 1/S.t.t; S.t.t := 1.0; S.t := S.t*D

  movsd   XMM8, One
  movhlps XMM9, XMM7    // XMM9 = S.t.t    |    w
  movlhps XMM7, XMM8    // XMM7 = S.t.z    |    S.t.t=1
  divsd   XMM8, XMM9    // XMM8 = 1.0/ S.z.z   = D
  movlhps XMM8, XMM8

  mulpd XMM6, XMM8
  mulpd XMM7, XMM8


  // save S.t
  movapd  oWORD[RCX+ 12*FloatSize], XMM6   //   S.z.x | S.z.y
  movapd  oWORD[RCX+ 14*FloatSize], XMM7   //   S.z.z |


// BACKWARD MOVEMENT!

// STEP  #7.x-y  D := S.x.y; S.x.y := 0; S.x -= S.y*D

  movhlps XMM8, XMM0     // XMM8= D=S.x.y | w
  movlhps XMM8, XMM8     // XMM8 = D | D

  movlhps XMM0, XMM10    // S.x.y = 0
  movapd XMM9, XMM2      // XMM9= S.y.x      | S.y.y
  mulpd XMM9, XMM8       // XMM9= S.y.x*D    | S.y.y*D
  subpd XMM0, XMM9       // XMM0= S.x.x+S.y.x*D | S.y.y*D

  movapd XMM9, XMM3       // XMM9= S.y.z      | t
  mulpd XMM9, XMM8       // XMM9= S.y.z*D    | t
  subpd XMM1, XMM9       // XMM0= S.x.z+S.y.z*D | t

// STEP  #7.y-z  D := S.y.z; S.y.z := 0; S.y -= S.z*D

  movsd   XMM8, XMM3     // XMM8= D=S.y.z | w
  movsd   XMM3, XMM10    // S.y.z := 0
  movlhps XMM8, XMM8     // XMM8 = D | D


  movapd XMM9, XMM4      // XMM9= S.z.x      | S.z.y
  mulpd XMM9, XMM8       // XMM9= S.z.x*D    | S.z.y*D
  subpd XMM2, XMM9       // XMM0= S.y.x-S.z.x*D | S.y.y-S.z.y*D

  movapd XMM9, XMM5       // XMM9= S.y.z      | t
  mulpd XMM9, XMM8       // XMM9= S.y.z*D    | t
  subpd XMM3, XMM9       // XMM0= S.x.z+S.y.z*D | t

// STEP  #7.x-z  D := S.x.z; S.x.z := 0; S.x -= S.z*D

  movsd   XMM8, XMM1     // XMM8= D=S.x.z | w
  movsd   XMM1, XMM10    // S.x.z := 0
  movlhps XMM8, XMM8     // XMM8 = D | D


  movapd XMM9, XMM4      // XMM9= S.t.x      | S.t.y
  mulpd XMM9, XMM8       // XMM9= S.t.x*D    | S.t.y*D
  subpd XMM0, XMM9       // XMM0= S.z.x-S.t.x*D | S.z.y-S.t.y*D

  movapd XMM9, XMM5      // XMM9= S.t.z      | S.t.t
  mulpd XMM9, XMM8       // XMM9= S.t.z*D    | S.t.t*D
  subpd XMM1, XMM9       // XMM0= S.z.z-S.t.z*D | S.z.t-S.t.t*D

// STEP  #7.z-t  D := S.z.t; S.z.t := 0; S.z -= S.t*D

  movhlps XMM8, XMM5     // XMM8= D=S.z.t | w
  movlhps XMM5, XMM10    // S.z.t := 0
  movlhps XMM8, XMM8     // XMM8 = D | D


  movapd XMM9, XMM6      // XMM9= S.t.x      | S.t.y
  mulpd XMM9, XMM8       // XMM9= S.t.x*D    | S.t.y*D
  subpd XMM4, XMM9       // XMM0= S.z.x-S.t.x*D | S.z.y-S.t.y*D

  movapd XMM9, XMM7      // XMM9= S.t.z      | S.t.t
  mulpd XMM9, XMM8       // XMM9= S.t.z*D    | S.t.t*D
  subpd XMM5, XMM9       // XMM0= S.z.z-S.t.z*D | S.z.t-S.t.t*D

  // SAVE S.z
  movupd oWORD[RCX+  8*FloatSize], XMM4   //   S.z.x | S.z.y
  movupd oWORD[RCX+ 10*FloatSize], XMM5   //   S.z.z |
// STEP  #8.y-t  D := S.y.t; S.y.t := 0; S.y -= S.t*D

  movhlps XMM8, XMM3     // XMM8= D=S.y.t | w
  movlhps XMM3, XMM10    // S.y.t = 0

  movlhps XMM8, XMM8     // XMM8 = D | D

  movapd XMM9, XMM6      // XMM9= S.y.x      | S.y.y
  mulpd XMM9, XMM8       // XMM9= S.y.x*D    | S.y.y*D
  subpd XMM2, XMM9       // XMM0= S.x.x+S.y.x*D | S.y.y*D

  movapd XMM9, XMM7      // XMM9= S.y.z      | t
  mulpd XMM9, XMM8       // XMM9= S.y.z*D    | t
  subpd XMM3, XMM9       // XMM0= S.x.z+S.y.z*D | t

//   save S.y
  movupd oWORD[RCX+  4*FloatSize], XMM2   //   S.y.x | S.y.y
  movupd oWORD[RCX+  6*FloatSize], XMM3   //   S.y.z |

// STEP  #8.x-t  D := S.x.t; S.x.t := 0; S.x -= S.t*D

  movhlps XMM8, XMM1     // XMM8= D=S.x.t | w
  movlhps XMM1, XMM10    // S.x.t := 0
  movlhps XMM8, XMM8     // XMM8 = D | D


  movapd XMM9, XMM6      // XMM9= S.t.x      | S.t.y
  mulpd XMM9, XMM8       // XMM9= S.t.x*D    | S.t.y*D
  subpd XMM0, XMM9       // XMM0= S.x.x-S.t.x*D | S.x.y-S.t.y*D

  movapd XMM9, XMM7      // XMM9= S.t.z         | S.t.t
  mulpd XMM9, XMM8       // XMM9= S.t.z*D       | S.t.t*D
  subpd XMM1, XMM9       // XMM0= S.x.z-S.t.z*D | S.x.t-S.t.t*D



  movupd oWORD[RCX              ], XMM0   //   S.x.x | S.x.y
  movupd oWORD[RCX+  2*FloatSize], XMM1   //   S.x.z |


  Add RCX, 16*FloatSize


  dec     EBX             // уменьшаем коунтер , значение падает в аккумулятор!
  jg      @next           // прыг на новую итерацию, если аккумулятор >0
{$ELSE}

var i: integer;
begin
  i := 0;
  repeat
    Invert_M4(  T_M4Arr(@S)[i] );
    inc(i);
  until i >= N;
{$ENDIF CPUX64}

end;


class function T_SSE.Invert_M3_clang(const M3: T_Tens3): T_Tens3;
const One : Uint64 = $3ff0000000000000 ;         //    # double 1
asm
//.LCPI0_0:
//        .quad   0x3ff0000000000000              # double 1
// Invert(Mat33):                        # @Invert(Mat33)
  .SAVENV XMM5
  .SAVENV XMM6
  .SAVENV XMM7
  .SAVENV XMM8
  .SAVENV XMM9
  .SAVENV XMM10
  .SAVENV XMM11
  .SAVENV XMM12

        mov     rax, rdi
        movsd   xmm8, qword ptr [M3 + 40-8]      //# xmm8 = mem[0],zero
        movsd   xmm12, qword ptr [M3 + 8-8]      //# xmm12 = mem[0],zero
        movsd   xmm9, qword ptr [M3 + 16-8]      //# xmm9 = mem[0],zero
        movsd   xmm4, qword ptr [M3 + 64-8]      //# xmm4 = mem[0],zero
        movapd  xmm10, xmm8
        unpcklpd        xmm10, xmm4                     //# xmm10 = xmm10[0],xmm4[0]
        unpcklpd        xmm4, xmm9                      //# xmm4 = xmm4[0],xmm9[0]
        movsd   xmm11, qword ptr [M3 + 72-8]            //# xmm11 = mem[0],zero
        movsd   xmm1, qword ptr [M3 + 48-8]             //# xmm1 = mem[0],zero
        movupd  xmm5, oWORD ptr [M3 + 24-8]     // movapd => movupd!
        movddup xmm7, xmm1                              //# xmm7 = xmm1[0,0]
        movapd  xmm13, xmm5
        unpcklpd        xmm13, xmm1                     //# xmm13 = xmm13[0],xmm1[0]
        unpcklpd        xmm1, xmm11                     //# xmm1 = xmm1[0],xmm11[0]
        mulpd   xmm1, xmm4
        movapd  xmm3, xmm8
        unpcklpd        xmm3, xmm11                     //# xmm3 = xmm3[0],xmm11[0]
        unpcklpd        xmm11, xmm5                     //# xmm11 = xmm11[0],xmm5[0]
        movapd  xmm6, xmm11
        mulpd   xmm6, xmm10
        movapd  xmm0, xmm6
        subpd   xmm0, xmm1
        movapd  xmm2, xmm0
        mulsd   xmm2, xmm12
        subpd   xmm1, xmm6
        movsd   xmm14, qword ptr [M3 + 32-8]     //# xmm14 = mem[0],zero
        mulpd   xmm1, xmm5
        unpckhpd        xmm1, xmm1              //# xmm1 = xmm1[1,1]
        subsd   xmm2, xmm1
        movsd   xmm1, qword ptr [M3 + 56-8]      //# xmm1 = mem[0],zero
        mulpd   xmm3, xmm5
        movapd  xmm6, xmm14
        mulsd   xmm14, xmm9
        unpcklpd        xmm9, xmm1                      //# xmm9 = xmm9[0],xmm1[0]
        mulpd   xmm7, xmm9
        subpd   xmm7, xmm3
        movapd  xmm3, xmm7
        mulsd   xmm3, xmm1
        addsd   xmm3, xmm2
        movsd   xmm2, One //qword ptr [rip + .LCPI0_0] //# xmm2 = mem[0],zero
        divsd   xmm2, xmm3
        movddup xmm2, xmm2                      //# xmm2 = xmm2[0,0]
        mulpd   xmm0, xmm2
        movupd  oWORD ptr [result], xmm0
        mulpd   xmm7, xmm2
        movupd  oWORD ptr [result + 16], xmm7
        unpcklpd        xmm6, xmm1                      //# xmm6 = xmm6[0],xmm1[0]
        unpcklpd        xmm1, xmm12                     //# xmm1 = xmm1[0],xmm12[0]
        mulpd   xmm13, xmm1
        blendpd xmm5, xmm12, 1                 // # xmm5 = xmm12[0],xmm5[1]
        mulpd   xmm5, xmm11
        subpd   xmm5, xmm13
        mulpd   xmm5, xmm2
        movupd  oWORD ptr [result + 32], xmm5
        mulpd   xmm1, xmm10
        mulpd   xmm6, xmm4
        subpd   xmm6, xmm1
        mulpd   xmm6, xmm2
        movupd  oWORD ptr [result + 48], xmm6
        mulsd   xmm12, xmm8
        subsd   xmm12, xmm14
        divsd   xmm12, xmm3
        movsd   qword ptr [result + 64], xmm12
end;

class function T_SSE.DotProd(const V0,V1 : T_Vect): real;
 {$IFDEF CPUX64}
asm
    .NOFRAME

    movapd XMM0, OWORD[V0]
    movapd XMM2, OWORD[V1]

    dppd XMM0, XMM2,  49 // imm8 = bit0 | ~bit1 | bit4 | bit5
                         // bit0- save res to low half
                         // bit4- * low parts
                         // bit5- * high parts

    movapd XMM1, OWORD[V0+2*FloatSize]
    movapd XMM3, OWORD[V1+2*FloatSize]

    mulsd XMM1, XMM3

    addsd XMM0, XMM1
 {$ELSE}
begin
  result := V0.x * V1.x + V0.y * V1.y+ V0.z * V1.z;
{$ENDIF CPUX64}


end;

class procedure T_SSE.DotProd(const V0,V1 : T_Vect; var Res: double);
 {$IFDEF CPUX64}
asm


    .NOFRAME
    movapd XMM0, OWORD[V0]
    movapd XMM2, OWORD[V1]

    dppd XMM0, XMM2,  49 // imm8 = bit0 | ~bit1 | bit4 | bit5
                         // bit0- save res to low half
                         // bit4- * low parts
                         // bit5- * high parts

    movapd XMM1, OWORD[V0+2*FloatSize]
    movapd XMM3, OWORD[V1+2*FloatSize]

    mulsd XMM1, XMM3

    addsd XMM0, XMM1
    movsd   QWORD[Res], xmm0


 {$ELSE}

begin
  Res := V0.x * V1.x + V0.y * V1.y+ V0.z * V1.z;
{$ENDIF CPUX64}

end;


class procedure T_SSE.Solve(const A: T_Tens; const B: T_Vect; var X: T_Vect);
{$IFDEF CPUX64}
const One : double = 1.0;
      Zero : double = 0.0;
asm

   .NOFRAME
   .SAVENV XMM4
   .SAVENV XMM5
   .SAVENV XMM6
   .SAVENV XMM7

   // Load matrix
   movapd  XMM0, oWORD[A]
   movapd  XMM1, oWORD[A+2 *FloatSize]
   movapd  XMM2, oWORD[A+4 *FloatSize]
   movapd  XMM3, oWORD[A+6 *FloatSize]
   movapd  XMM4, oWORD[A+8 *FloatSize]
   movapd  XMM5, oWORD[A+10*FloatSize]

   // Load FreePart to extended row of matrix
   movhpd  XMM1, qWORD[B             ]
   movhpd  XMM3, qWORD[B+ FloatSize  ]
   movhpd  XMM5, qWORD[B+ 2*FloatSize]


   // divide first string with a00 ( make 1.0 on main diagonal)
{   movlpd  XMM7, One
   divsd   XMM7, XMM0    // XMM7 = 1/A00 | 0
   movlhps XMM7, XMM7    // XMM7 = 1/A00 | 1/A00

   mulpd   XMM0, XMM7
   mulpd   XMM1, XMM7    // XMM0, XMM1 = 1.0 | A01/A00,  A02/A00 | b0/A00 }


   movapd  XMM7, XMM0
   movlhps XMM7, XMM7

   movlpd  XMM0, One

   divpd   XMM0, XMM7    // XMM0 = 1/A00 | A01/A00
   movapd  XMM7, XMM0    // XMM7 = 1/A00 | 1/A00
   movlhps XMM7, XMM0

   mulpd   XMM1, XMM7    // XMM0, XMM1 = 1.0/a00 | a01/a00,  a02/a00 | b0/a00


   // get elem of first string, multiply with a10 and subtract from second string
   movsd   XMM7, XMM2
   movlhps XMM7, XMM7    // XMM7 = a10 | a10

   movapd  XMM6, XMM0    // use XMM6 as temp accumulator
   mulpd   XMM6, XMM7
   subpd   XMM2, XMM6

   movapd  XMM6, XMM1
   mulpd   XMM6, XMM7
   subpd   XMM3, XMM6    // XMM2, XMM3 = 0 | A11*,  A12* | b1*

   // get elem of first string, multiply with a20 and subtract from third string
   movsd   XMM7, XMM4
   movlhps XMM7, XMM7    // XMM7 = A20 | A20

   movapd  XMM6, XMM0
   mulpd   XMM6, XMM7
   subpd   XMM4, XMM6

   movapd  XMM6, XMM1
   mulpd   XMM6, XMM7
   subpd   XMM5, XMM6    // XMM4, XMM5 = 0 | A21*,  A22* | b2*

   movhlps XMM2,XMM2     // XMM2 = a11* | a11*
   divpd XMM3, XMM2      // XMM3 = a12** | b1**

   movhlps XMM4, XMM4    // XMM4 = a21 | a21
   mulpd   XMM4, XMM3
   subpd   XMM5, XMM4    // XMM5 = a22*, b2*
   movhlps XMM7, XMM5    // XMM7 = b2*
   divsd   XMM7, XMM5    // XMM7 = b2*/ a22*   = x2

   movapd  XMM4, XMM1    // XMM4 = a02 | b0
   movlhps XMM4, XMM3    // XMM4 = a02 | a12

   movlhps XMM7, XMM7    // XMM7 = x2 | x2

   mulpd   XMM4, XMM7    // XMM4 = a02*x2 | a12*x2
   movhlps XMM3, XMM1    // XMM3 = b0* | b1*

   subpd   XMM3, XMM4    // XMM3 = b0** | x1
   mulpd   XMM0, XMM3
   movhlps XMM0, XMM0
   subsd   XMM3, XMM0

   movapd  oWORD[X              ],  XMM3
   movsd   qWORD[X + 2*FloatSize],  XMM7

{$ELSE}
var D: real;
    r: T_Vect;
begin
  D:=    A.x.x*( A.y.y*A.z.z - A.y.z*A.z.y )
       + A.x.y*( A.z.x*A.y.z - A.y.x*A.z.z )
       + A.x.z*( A.y.x*A.z.y - A.y.y*A.z.x ) ;

  if Abs(D)> 1.0e-16 then
    D := 1.0 / D
  else begin
    X :=  B;
    exit;
  end;

  r := B;

  X.x :=(    (  A.y.y*A.z.z - A.y.z*A.z.y  ) * r.x
            -(  A.x.y*A.z.z - A.x.z*A.z.y  ) * r.y
            +(  A.x.y*A.y.z - A.x.z*A.y.y  ) * r.z  ) * D;

  X.y :=(   -(  A.y.x*A.z.z - A.z.x*A.y.z  ) * r.x
            +(  A.x.x*A.z.z - A.x.z*A.z.x  ) * r.y
            -(  A.x.x*A.y.z - A.y.x*A.x.z  ) * r.z  ) * D;

  X.z :=(    (  A.y.x*A.z.y - A.y.y*A.z.x  ) * r.x
            -(  A.x.x*A.z.y - A.x.y*A.z.x  ) * r.y
            +(  A.x.x*A.y.y - A.x.y*A.y.x  ) * r.z  ) * D;
{$ENDIF}

end;


class procedure T_SSE.Solve(const A: P_Tens; const B: P_Vect; const X: P_Vect;
    const N: integer);
{$IFDEF CPUX64}
const One : double = 1.0;
      Zero : double = 0.0;
asm

   .NOFRAME
   .SAVENV XMM4
   .SAVENV XMM5
   .SAVENV XMM6
   .SAVENV XMM7

   .PUSHNV RBX

   mov EBX, N
   // Load matrix

 @next:

   movapd  XMM0, oWORD[A             ]
   movapd  XMM1, oWORD[A+2 *FloatSize]
   movapd  XMM2, oWORD[A+4 *FloatSize]
   movapd  XMM3, oWORD[A+6 *FloatSize]
   movapd  XMM4, oWORD[A+8 *FloatSize]
   movapd  XMM5, oWORD[A+10*FloatSize]

   // Load FreePart to extended row of matrix
   movhpd  XMM1, qWORD[B             ]
   movhpd  XMM3, qWORD[B+ FloatSize  ]
   movhpd  XMM5, qWORD[B+ 2*FloatSize]


   // divide first string with a00 ( make 1.0 on main diagonal)
{   movlpd  XMM7, One
   divsd   XMM7, XMM0    // XMM7 = 1/A00 | 0
   movlhps XMM7, XMM7    // XMM7 = 1/A00 | 1/A00

   mulpd   XMM0, XMM7
   mulpd   XMM1, XMM7    // XMM0, XMM1 = 1.0 | A01/A00,  A02/A00 | b0/A00 }


   movapd  XMM7, XMM0   // XMM7 = a00 | trash
   movlhps XMM7, XMM7   // XMM7 = a00 | a00

   movlpd  XMM0, One     // XMM0 = 1.0 | a01

   divpd   XMM0, XMM7    // XMM0 = 1/A00 | A01/A00
   movapd  XMM7, XMM0    // XMM7 = 1/A00 | 1/A00
   movlhps XMM7, XMM0

   mulpd   XMM1, XMM7    // XMM0, XMM1 = 1.0/a00 | a01/a00,  a02/a00 | b0/a00


   // get elem of first string, multiply with a10 and subtract from second string
   movsd   XMM7, XMM2    // XMM7 = a10 | trash
   movlhps XMM7, XMM7    // XMM7 = a10 | a10

   movapd  XMM6, XMM0    // use XMM6 as temp accumulator
   mulpd   XMM6, XMM7
   subpd   XMM2, XMM6

   movapd  XMM6, XMM1
   mulpd   XMM6, XMM7
   subpd   XMM3, XMM6    // XMM2, XMM3 = 0 | A11*,  A12* | b1*

   // get elem of first string, multiply with a20 and subtract from third string
   movsd   XMM7, XMM4
   movlhps XMM7, XMM7    // XMM7 = A20 | A20

   movapd  XMM6, XMM0
   mulpd   XMM6, XMM7
   subpd   XMM4, XMM6

   movapd  XMM6, XMM1
   mulpd   XMM6, XMM7
   subpd   XMM5, XMM6    // XMM4, XMM5 = 0 | A21*,  A22* | b2*

   movhlps XMM2, XMM2    // XMM2 = a11* | a11*
   divpd   XMM3, XMM2    // XMM3 = a12** | b1**

   movhlps XMM4, XMM4    // XMM4 = a21 | a21
   mulpd   XMM4, XMM3
   subpd   XMM5, XMM4    // XMM5 = a22*, b2*
   movhlps XMM7, XMM5    // XMM7 = b2*
   divsd   XMM7, XMM5    // XMM7 = b2*/ a22*   = x2

   movapd  XMM4, XMM1    // XMM4 = a02 | b0
   movlhps XMM4, XMM3    // XMM4 = a02 | a12

   movlhps XMM7, XMM7    // XMM7 = x2 | x2

   mulpd   XMM4, XMM7    // XMM4 = a02*x2 | a12*x2
   movhlps XMM3, XMM1    // XMM3 = b0* | b1*

   subpd   XMM3, XMM4    // XMM3 = b0** | x1
   mulpd   XMM0, XMM3
   movhlps XMM0, XMM0
   subsd   XMM3, XMM0

   movapd  oWORD[X              ],  XMM3
   movsd   qWORD[X + 2*FloatSize],  XMM7

   Add A, 12*FloatSize
   Add B, 4*FloatSize
   Add X, 4*FloatSize

   dec     EBX             // уменьшаем коунтер , значение падает в аккумулятор!
   jg     @next          // прыг на новую итерацию, если аккумулятор <>0



{$ELSE}
var i: integer;
    AA: P_Tens;
    BB: P_Vect;
    XX: P_Vect;
begin
  AA := A;
  BB := B;
  XX := X;
  for i := 0 to N-1 do
  begin
//    Solve_old( AA^, BB^, XX^ );
    Gauss( AA^, BB^, XX^ );
    Inc(AA);
    Inc(BB);
    Inc(XX);// := XX + 48;

  end;
{$ENDIF}

end;


class procedure T_SSE.Solve_old(const A: T_Tens; const B: T_Vect; var X:
    T_Vect);
{$IFDEF CPUX64}
const One : double = 1.0;
      Zero : double = 0.0;
asm

// approximately 10% slower the best current result
// ~177 ms on 768MB array at Core-i5-4460.
// (fastest result ~145 ms on the same pc.

   .SAVENV XMM4
   .SAVENV XMM5
   .SAVENV XMM6
   .SAVENV XMM7
   movapd XMM0, oWORD[A]
   movapd XMM1, oWORD[A+2 *FloatSize]
   movapd XMM2, oWORD[A+4 *FloatSize]
   movapd XMM3, oWORD[A+6 *FloatSize]
   movapd XMM4, oWORD[A+8 *FloatSize]
   movapd XMM5, oWORD[A+10*FloatSize]


   movhpd XMM1, qWORD[B             ]
   movhpd XMM3, qWORD[B+ FloatSize  ]
   movhpd XMM5, qWORD[B+ 2*FloatSize]

   movlpd XMM7, One
   divsd XMM7, XMM0    // XMM7 = 1/A00 | 0
   movlhps XMM7, XMM7  // XMM7 = 1/A00 | 1/A00

   mulpd XMM0, XMM7    // XMM0 = 1.0, A01/A00, A02/A00, b0/A00
   mulpd XMM1, XMM7

   movsd XMM7, XMM2
   movlhps XMM7, XMM2  // XMM7 = A10 | A10

   movapd XMM6, XMM0
   mulpd XMM6, XMM7
   subpd XMM2, XMM6

   movapd XMM6, XMM1
   mulpd XMM6, XMM7
   subpd XMM3, XMM6


   movsd XMM7, XMM4
   movlhps XMM7, XMM4  // XMM7 = A20 | A20

   movapd XMM6, XMM0
   mulpd XMM6, XMM7
   subpd XMM4, XMM6

   movapd XMM6, XMM1
   mulpd XMM6, XMM7
   subpd XMM5, XMM6


   movlpd  XMM7, One  // XMM7 = 1.0 | trahs
   movhlps XMM6, XMM2 // XMM6 = a11* | trash
   divsd   XMM7, XMM6 // XMM7 = 1.0/a11
   movlhps XMM7, XMM7 // XMM7 = 1/a11* | 1/a11*

   mulpd XMM3, XMM7  // XMM3 = a12/a11, b1/a11

   movhlps XMM7, XMM4  // XMM7 = a21
   movlhps XMM7, XMM7  // XMM7 = a21 | a21

   mulpd XMM7, XMM3
   subpd XMM5, XMM7    // XMM5 = a22*, b2*
   movhlps XMM7, XMM5  // XMM7 = b2*
   divsd XMM7, XMM5    // XMM7 = x2 = b2*/ a22*

   mulsd XMM3, XMM7    // XMM3 = a12*x2 | b1*
   mulsd XMM1, XMM7    // XMM1 = a02*x2 | b0*

   movhlps XMM6, XMM3   // XMM6= b1*
   subsd XMM6, XMM3     // XMM6 = b1* - a12*x2
   movhlps XMM0, XMM0   // XMM0 = a01*

   mulsd XMM0, XMM6     // XMM0 = a01*x1

   addsd XMM0, XMM1     // XMM0 = a01*x1 + a02*x2
   movhlps XMM1, XMM1   // XMM1 = b0* | b0*

   subsd XMM1, XMM0     // XMM1 = b0* - a01*x1 - a02*x2 = x0!

   movsd qWORD[X              ], XMM1
   movsd qWORD[X +   FloatSize], XMM6
   movsd qWORD[X + 2*FloatSize], XMM7

{$ELSE}
var D: real;
    r: T_Vect;
begin
  D:=    A.x.x*( A.y.y*A.z.z - A.y.z*A.z.y )
       + A.x.y*( A.z.x*A.y.z - A.y.x*A.z.z )
       + A.x.z*( A.y.x*A.z.y - A.y.y*A.z.x ) ;

  if Abs(D)> 1.0e-16 then
    D := 1.0 / D
  else begin
    X :=  B;
    exit;
  end;

  r := B;

  X.x :=(    (  A.y.y*A.z.z - A.y.z*A.z.y  ) * r.x
            -(  A.x.y*A.z.z - A.x.z*A.z.y  ) * r.y
            +(  A.x.y*A.y.z - A.x.z*A.y.y  ) * r.z  ) * D;

  X.y :=(   -(  A.y.x*A.z.z - A.z.x*A.y.z  ) * r.x
            +(  A.x.x*A.z.z - A.x.z*A.z.x  ) * r.y
            -(  A.x.x*A.y.z - A.y.x*A.x.z  ) * r.z  ) * D;

  X.z :=(    (  A.y.x*A.z.y - A.y.y*A.z.x  ) * r.x
            -(  A.x.x*A.z.y - A.x.y*A.z.x  ) * r.y
            +(  A.x.x*A.y.y - A.x.y*A.y.x  ) * r.z  ) * D;
{$ENDIF}

end;


class procedure T_SSE.DotProd(const V0,V1 : T_VectArr;const Res: T_RealArr; const
    i0,i1: integer);
 {$IFDEF CPUX64}
asm

    .NOFRAME

    XOR RAX, RAX
    mov EAX, i0;

    imul RAX, FloatSize;

    add Res, RAX
    imul RAX, 4
    add V0, RAX
    add V1, RAX

    XOR RAX, RAX
    mov EAX, i1
    sub EAX, i0
    inc RAX

 @next:
    movapd XMM0, OWORD[V0]
    movapd XMM2, OWORD[V1]

    dppd XMM0, XMM2,  49 // imm8 = bit0 | ~bit1 | bit4 | bit5
                         // bit0- save res to low half
                         // bit4- * low parts
                         // bit5- * high parts

    movapd XMM1, OWORD[V0+2*FloatSize]
    movapd XMM3, OWORD[V1+2*FloatSize]

    mulsd XMM1, XMM3

    addsd XMM0, XMM1
    movsd   QWORD[Res], xmm0

    add V0, VectSize
    add V1, VectSize
    add Res, FloatSize

  dec     RAX             // уменьшаем коунтер , значение падает в аккумулятор!
  jg      @next           // прыг на новую итерацию, если аккумулятор >0

{$ELSE}
begin

{$ENDIF CPUX64}

end;







class procedure T_SSE.DotProd(const V0,V1 : T_VectArr; var Res: real; const i0,
    i1: integer);
 {$IFDEF CPUX64}
asm

    .NOFRAME

    XOR RAX, RAX
    mov EAX, i0;

    imul RAX, FloatSize;

    add Res, RAX
    imul RAX, 4
    add V0, RAX
    add V1, RAX

    XOR RAX, RAX
    mov EAX, i1
    sub EAX, i0
    inc RAX

    PXOR XMM3, XMM3
 @next:
    movapd XMM0, OWORD[V0]
    movapd XMM1, OWORD[V1]

    mulpd  XMM0, XMM1

    addpd XMM3, XMM0



    movapd XMM0, OWORD[V0+2*FloatSize]
    movapd XMM1, OWORD[V1+2*FloatSize]

    mulsd XMM0, XMM1


    addsd XMM3, XMM0

    add V0, VectSize
    add V1, VectSize


  dec     RAX             // уменьшаем коунтер , значение падает в аккумулятор!
  jg      @next           // прыг на новую итерацию, если аккумулятор >0

  movhlps XMM0, XMM3
  addsd XMM0, XMM3

  movsd qWORD[Res], XMM0


{$ELSE}
begin

{$ENDIF CPUX64}

end;




function T_Tens3.ToTens: T_Tens;
begin
  Result.x.x := xx;
  Result.x.y := xy;
  Result.x.z := xz;


  Result.y.x := yx;
  Result.y.y := yy;
  Result.y.z := yz;

  Result.z.x := zx;
  Result.z.y := zy;
  Result.z.z := zz;


end;

class operator T_Tens3.Explicit(A: T_Tens3): T_Tens;
begin
  Result.x.x := A.xx;
  Result.x.y := A.xy;
  Result.x.z := A.xz;


  Result.y.x := A.yx;
  Result.y.y := A.yy;
  Result.y.z := A.yz;

  Result.z.x := A.zx;
  Result.z.y := A.zy;
  Result.z.z := A.zz;


end;

class operator T_Tens3.Implicit(A: T_Tens): T_Tens3;
begin
  Result.xx := A.x.x;
  Result.xy := A.x.y;
  Result.xz := A.x.z;


  Result.yx := A.y.x;
  Result.yy := A.y.y;
  Result.yz := A.y.z;

  Result.zx := A.z.x;
  Result.zy := A.z.y;
  Result.zz := A.z.z;


end;


end.


