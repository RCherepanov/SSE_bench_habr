unit sse_Main;

interface

(*
Тестовый проект для сравнения скорости работы различных реализаций векторной алгебры
Базовый набор функций лежат в base.types.pas

Для проверки скорости работы какой-либо функции нужно сделать следующее:

for i := 0 to N-1 do
  Заполнить_значениями_нужные_массивы(X, Xsse, V1,V2,V3, M3x3, M3x3compact, M4x4)

Start;
выполнить_нужные_операции_над_массивами
Stop;
WriteResult( 'FunctionName', SizeOfReadDataPerIteration );

SizeOfReadDataPerIteration сумма размеров элементов всех массивов, которые считывались при обработке.
то есть, если при обработке читались V2,V3, M3x3, а результат записывался в V2- то
SizeOfReadDataPerIteration = SizeOf(V2[0])+ SizeOf(V3)  + SizeOf(M3x3)
 если при обработке читались V3, M3x3, а результат записывался в V2- то V2 не учитываем, и
SizeOfReadDataPerIteration = SizeOf(V3)  + SizeOf(M3x3)

если при выводе результатов хочется убедиться в корректности какого-либо значения- то можно использовать
переопределенные варианты WriteResult


При вызове метода Start вызывается preheat_CPU- некоторое количество расчетов
логарифма- симулирует нагрузку для того, чтобы вывести CPU на повышенную частоту.
без этого результаты сильно скачут.

Черепанов Р.О.     08 March 2021
RCherepanov82@gmail.com
Tiriet@habr.com

*)


uses
  FastMM4,
  base.types,
  miniprof,
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.Menus;

type
  TForm1 = class(TForm)
    bnt_Mij_Vj: TButton;
    bnt_ViBij: TButton;
    btn_DotProd_Array: TButton;
    btn_DotProd_summ: TButton;
    btn_TestFastInvert: TButton;
    btn_TestFastInvert_4x4: TButton;
    est1: TMenuItem;
    File1: TMenuItem;
    MainMenu1: TMainMenu;
    Memo1: TMemo;
    procedure bnt_Mij_VjClick(Sender: TObject);
    procedure bnt_ViBijClick(Sender: TObject);
    procedure btn_DotProd_ArrayClick(Sender: TObject);
    procedure btn_DotProd_summClick(Sender: TObject);
    procedure btn_TestFastInvertClick(Sender: TObject);
    procedure btn_TestFastInvert_4x4Click(Sender: TObject);
    procedure FormCreate(Sender: TObject);
  private
    ArraySize: Int64;
    CalcTime: Int64;
    fCalcTimeToSec: Double;
    fStartTime: Int64;
    fStopTime: Int64;
    M3x3: T_TensorArr;
    M3x3compact: T_TensorArr3;
    M4x4: T_M4Arr;
    V1: T_VectArr;
    V2: T_VectArr;
    V3: T_VectArr;
    X: T_RealArr;
    Xsse: T_RealArr;
    function GetCalcTime: Int64;
    procedure WriteResult(const expected, actual: real; const FuncName: string;
        const DataSize: integer); overload;
    procedure WriteResult(const FuncName: string; const DataSize: integer);
        overload;
    procedure WriteResult(const base, inverted: T_M4; const FuncName: string; const
        DataSize: integer); overload;
    procedure WriteResult(const base, inverted: T_Tens; const FuncName: string;
        const DataSize: integer); overload;
    procedure _Write(const S: string; const Params: array of const);
    procedure _WriteMatrix(const header: String; const m: T_M4); overload;
    procedure _WriteMatrix(const header: String; const m: T_Tens); overload;
    procedure _WriteVector(const header: String; const V3: T_Vect); overload;
  public
    procedure preheat_cpu;
    procedure Start;
    procedure Stop;
  published
    procedure Test_DotProd_arrays;
    procedure Test_DotProd_summ;
    procedure Test_Inversion_3x3;
    procedure Test_Mij_Vj;
    procedure Test_Vj_Mji;
  end;

var
  Form1: TForm1;

implementation

{$R *.dfm}
const N: integer = 8*1024*1024;


procedure TForm1.bnt_Mij_VjClick(Sender: TObject);
begin
  Test_Mij_Vj;
end;

procedure TForm1.bnt_ViBijClick(Sender: TObject);
begin
  Test_Vj_Mji;
end;

procedure TForm1.btn_DotProd_ArrayClick(Sender: TObject);
begin
  Test_DotProd_arrays;
end;

procedure TForm1.btn_DotProd_summClick(Sender: TObject);
begin
  Test_DotProd_summ;
end;

procedure TForm1.btn_TestFastInvertClick(Sender: TObject);
begin
  Test_Inversion_3x3;
end;

procedure TForm1.btn_TestFastInvert_4x4Click(Sender: TObject);
var t1, t2: T_M4;
    t_0, t2orig: T_M4;
    err_orig, err_fast: real;
    i: integer;




begin


  Memo1.Lines.Add(Format('1s = %d calcTime',[CalcTime]));

//  randseed := 0;
  t_0 := Matrix4(
     random +1.0, random -0.5, random -0.5, random -0.5,
     random -0.5, random +1.0, random -0.5, random -0.5,
     random -0.5, random -0.5, random +1.0, random -0.5,
     random -0.5, random -0.5, random -0.5, random +1.0
  );

  t1 := t_0;
  t2 := t_0;

  Invert(t1);

  _WriteMatrix( 'T^-1 (orig) = ', t1);
  Memo1.Lines.Add( '' );
  t1 := mult(t1, t_0);


  if Norma(t1, M4_Unite) > 1.0e-14 then
  begin
      _Write( 'residual',[] );
      _WriteMatrix( 'T*T^-1 (orig) = ', t1);
  end
  else
    _Write('Residual: %4.3g', [Norma(t1, M4_unite)]);
  Memo1.Lines.Add( '' );

  t1 := t_0;

  T_SSE.Invert_gauss(t1);

  _WriteMatrix( 'T^-1 (sse.old) = ', t1);
  Memo1.Lines.Add( '' );
  t1 := mult(t1, t_0);

  if Norma(t1, M4_Unite) > 1.0e-14 then
  begin
    _Write( 'residual',[] );
    _WriteMatrix( 'T*T^-1 (orig) = ', t1);
  end
  else
    _Write('Residual: %4.3g', [Norma(t1, M4_unite)]);
  Memo1.Lines.Add( '' );

/////////////////////

  Memo1.Lines.Add('');

/////////////////////////
///
  for i := 0 to N-1 do    M4x4[i] := t_0;

  Start;
  for i := 0 to N-1 do    Invert(M4x4[i]);
  Stop;

  WriteResult(t_0,  M4x4[2], 'M4: Invert.original', SizeOf(t_0));

/////////////////////////
  for i := 0 to N-1 do    M4x4[i] := t_0;

  Start;
  for i := 0 to N-1 do    T_SSE.Invert_gauss(M4x4[i]);
  Stop;

  WriteResult(t_0,  M4x4[2], 'M4: Invert_SSE_gauss', SizeOf(t_0));

/////////////////////////
  for i := 0 to N-1 do    M4x4[i] := t_0;

  Start;
    T_SSE.Invert_gauss(M4x4[0], N);
  Stop;

  WriteResult(t_0,  M4x4[2], 'M4: Invert_SSE_gauss(T_M4, N)', SizeOf(t_0));

  SetLength(M4x4, 0);
end;

procedure TForm1.FormCreate(Sender: TObject);
var i: integer;
    d1,d2: pointer;
    t1,t2: TDateTime;
begin

  t1 := Now;
  start;

  sleep(100);

  t2 := Now;
  stop;

  t2 := t2-t1;

  {$IFDEF CPUX64}
    Caption := Caption +' (X64)';
  {$ELSE}
    Caption := Caption +' (X32)';
  {$ENDIF}


  fCalcTimeToSec := t2*24*60*60/CalcTime;


  SetMinimumBlockAlignment(mba16Byte);

  GetMem(d1, 2);
  GetMem(d2, 2);

  SetLength(V3, N+1);
  SetLength(V1, N+1);
  SetLength(V2, N+1);
  SetLength(M3x3, N+1);
  SetLength(M3x3compact, N+1);
  SetLength(M4x4, N+1);
  SetLength(X, N+1);
  SetLength(Xsse, N+1);


  RANDSEED:= (0);
  for i := 0 to N-1 do
  begin
    V3[i] := Vect(Random, random, random);
    V1[i] := Vect(Random, random, random);
    V2[i] := Vect(Random, random, random);

    M3x3[i] := Tensor(
      1.0, 0.3, 0.2,
      0.3, 1.0, 0.1,
      0.2, 0.1, 1.0
    );

  end;

  Memo1.Lines.Clear;
  Memo1.Lines.Add(   Format( 'Array Aligment %d', [integer(@(V3[0])) mod 16]  ));
  Memo1.Lines.Add(   Format( 'Pointer Aligment %d, %d',
      [integer(d1) mod 16, integer(d2) mod 16]  ));

  FreeMem(d1);
  FreeMem(d2);

end;

function TForm1.GetCalcTime: Int64;
var Freq: Int64;
begin

  QueryPerformanceFrequency(Freq);
  Result := ((fStopTime- fStartTime)*1000000) div Freq ;

end;

procedure TForm1.WriteResult(const expected, actual: real; const FuncName: string;
    const DataSize: integer);
begin

   if Abs( expected - actual )>  1.0e-8*(1.0+Abs(Expected)) then begin
     _Write('%S: ERROR, expected: %4.3f; actual %4.3f', [FuncName, Expected, Actual]);
     exit;
   end;

   WriteResult( FuncName, DataSize ) ;
end;

procedure TForm1.WriteResult(const FuncName: string; const DataSize: integer);
begin
   _Write('"%S": CalcTime =%d mks; throutput %4.1f MB/s',
         [FuncName, CalcTime, ArraySize*DataSize / CalcTime ] ) ;
end;

procedure TForm1.WriteResult(const base, inverted: T_M4; const FuncName: string;
    const DataSize: integer);
    var Test: T_M4;
begin
   Test := Mult(base, Inverted);
   Add(Test, -1.0);

   if Norma(Test)> 1.0e-10 then begin
     _Write('%S: ERROR, Base*Inverted <> I', [FuncName]);
     exit;
   end;

   WriteResult( FuncName, DataSize ) ;
end;

procedure TForm1.WriteResult(const base, inverted: T_Tens; const FuncName:
    string; const DataSize: integer);
var Test: T_Tens;
begin
   Test := Mult(base, Inverted);
   Add(Test, -1.0);

   if Norma(Test)> 1.0e-10 then begin
     _Write('%S: ERROR, Base*Inverted <> I', [FuncName]);
     exit;
   end;

   WriteResult( FuncName, DataSize ) ;
end;

procedure TForm1._Write(const S: string; const Params: array of const);
begin
  Memo1.Lines.Add( Format(S, Params));
end;

procedure TForm1._WriteMatrix(const header: String; const m: T_M4);
const FmtStr : string = '%9.8g; %9.8g; %9.8g;  %9.8g; ';
begin
  Memo1.Lines.Add( header );
  Memo1.Lines.Add( Format( FmtStr, [m[0][0], m[0][1], m[0][2], m[0][3]] )   );
  Memo1.Lines.Add( Format( FmtStr, [m[1][0], m[1][1], m[1][2], m[1][3]] )   );
  Memo1.Lines.Add( Format( FmtStr, [m[2][0], m[2][1], m[2][2], m[2][3]] )   );
  Memo1.Lines.Add( Format( FmtStr, [m[3][0], m[3][1], m[3][2], m[3][3]] )   );
end;

    procedure TForm1._WriteMatrix(const header: String; const m: T_Tens);
    const FmtStr : string = '%9.8g; %9.8g; %9.8g; ';
    begin
      Memo1.Lines.Add( header );
      Memo1.Lines.Add( Format( FmtStr, [m.x.x, m.x.y, m.x.z] )   );
      Memo1.Lines.Add( Format( FmtStr, [m.y.x, m.y.y, m.y.z] )   );
      Memo1.Lines.Add( Format( FmtStr, [m.z.x, m.z.y, m.z.z] )   );
    end;

    procedure TForm1._WriteVector(const header: String; const V3: T_Vect);
    const FmtStr : string = '%9.8g; %9.8g; %9.8g; ';
    begin
      Memo1.Lines.Add( header + '= '+ Format( FmtStr, [V3.x, V3.y, V3.z] )   );
    end;

procedure TForm1.preheat_cpu;
var i: integer;
    s: real;
begin
  // some hard calculation to Boost processor clock
  s := 0.0;
  for i := 0 to 1024*1024 do
    s := s + ln(1.0+i/1000);

  if s<0.5 then MessageDLG('ln sum error', mtError, [mbOk],0);

end;

procedure TForm1.Start;
begin
  ArraySize := N;


  preheat_cpu;
  QueryPerformanceCounter(fStartTime);

end;

procedure TForm1.Stop;
begin
  QueryPerformanceCounter(fStopTime);
  CalcTime := GetCalcTime;
end;

procedure TForm1.Test_DotProd_arrays;
var i: integer;
    Res, Res_sse: real;

begin
  Memo1.Lines.Add('');

  Start;
  Res := 0;
  for i := 0 to N-1 do
    X[i] := DotProd(V1[i], V2[i]);
  Stop;

  for i := 0 to N-1 do
    Res := Res + X[i];

  WriteResult( Res, Res, 'Array of dot.prod: V1[i]*V2[i]', 2*SizeOf(V1[0]) );


  for i := 0 to N-1 do
    X[i] := 0.0;


  Start;
  for i := 0 to N-1 do
    X[i] := T_SSE.DotProd(V1[i], V2[i]);
  Stop;

  Res_SSE := 0;
  for i := 0 to N-1 do
    Res_SSE := Res_SSE + X[i];

  WriteResult( Res, Res_SSE, 'Array of dot.prod SSE: V1[i]*V2[i]', 2*SizeOf(V3[0]) );


  for i := 0 to N-1 do
    X[i] := 0.0;


  Start;
    T_SSE.DotProd(V1, V2, X, 0, N-1 );
  Stop;

  Res_SSE := 0;
  for i := 0 to N-1 do
    Res_SSE := Res_SSE + X[i];

  WriteResult( Res, Res_SSE, 'Array of dot.prod SSE, loop: V1[i]*V2[i]', 2*SizeOf(V3[0]) );

end;

procedure TForm1.Test_DotProd_summ;
var i: integer;
    Res, Res_sse: real;

begin
  Memo1.Lines.Add('');

  SetLength(X, N);

  Res := 0;

  Start;
  for i := 0 to N-2 do    Res := Res + DotProd(V3[i], V3[i+1]);
  Stop;

  WriteResult( res, res, 'Summ of dot.prod compiler V3[i]*V3[i+1]', 2*SizeOf(V3[0]) );


  Res_SSE := 0;

  Start;
  for i := 0 to N-2 do    Res_SSE := Res_SSE + T_SSE.DotProd(V3[i], V3[i+1]);
  Stop;

  WriteResult( res, Res_SSE, 'Summ of dot.prod SSE V3[i]*V3[i+1]', 2*SizeOf(V3[0]) );



  Start;
  for i := 0 to N-1 do    Res := Res + DotProd(V1[i], V2[i]);
  Stop;

  WriteResult( res, res, 'Summ of dot.prod compiler V1[i]*V2[i]', 2*SizeOf(V3[0]) );



  Start;
  Res_SSE := 0;
  for i := 0 to N-1 do
    Res_SSE := Res_SSE + T_SSE.DotProd(V1[i], V2[i]);
  Stop;

  WriteResult( Res_SSE, Res_SSE, 'Summ of dot.prod SSE: V1[i]*V2[i]', 2*SizeOf(V1[0]) );
  Res := Res_SSE;


  Start;
  Res_SSE := 0;
    T_SSE.DotProd(V1, V2, Res_SSE, 0, N-1 );
  Stop;

  WriteResult( Res, Res_SSE, 'Summ of dot.prod SSE, loop: V1[i]*V2[i]', 2*SizeOf(V1[0]) );
  Res := Res_SSE;



end;

procedure TForm1.Test_Inversion_3x3;
var t1, t2: T_Tens;
    t_0, t2orig: T_Tens;
    err_orig, err_fast: real;
    i: integer;

begin
//  randseed := 0;
  t_0 := Tensor(
     random +1.0, random -0.5, random +2,
     random -0.5, random +1.0, random -0.5,
     random +2.5, random -0.5, random +1.0
  );

  t_0 := Tensor(
     1, 2, 3,
     5, 6, 7,
     4, 9, 8
  );

  t1 := t_0;
  t2 := t_0;


  Invert (t1);

  _WriteMatrix( 'T^-1 (orig) = ', t1);
  Memo1.Lines.Add( '' );
  t1 := t1*t_0;


  if Norma(t1, S1) > 1.0e-14 then
  begin
      _Write( 'residual',[] );
      _WriteMatrix( 'T*T^-1 (orig) = ', t1);
  end
  else
    _Write('Residual: %4.3g', [Norma(t1, S1)]);
  Memo1.Lines.Add( '' );

  t1 := t_0;

  Invert_gauss(t1);

  _WriteMatrix( 'T^-1 (fast) = ', t1);
  Memo1.Lines.Add( '' );
  t1 := t1*t_0;

  if Norma(t1, S1) > 1.0e-14 then
  begin
      _Write( 'residual',[] );
      _WriteMatrix( 'T*T^-1 (orig) = ', t1);
  end
  else
    _Write('Residual: %4.3g', [Norma(t1, S1)]);
  Memo1.Lines.Add( '' );

  t1 := t_0;
{  t1 := Tensor(
     1, 2, 3,
     1, 2, 3,
     4, 9, 8
  ); }

  T_SSE.Invert(t1);

  _WriteMatrix( 'T^-1 (SSE) = ', t1);
  Memo1.Lines.Add( '' );
  t1 := t1*t_0;
  if Norma(t1, S1) > 1.0e-14 then
  begin
      _Write( 'residual',[] );
      _WriteMatrix( 'T*T^-1 (orig) = ', t1);
  end
  else
    _Write('Residual: %4.3g', [Norma(t1, S1)]);
  Memo1.Lines.Add( '' );


/////////////////////

  Memo1.Lines.Add('');


  for i := 0 to N-1 do    M3x3[i] := t_0;
  Start;
  for i := 0 to N-1 do
    Invert(M3x3[i]);
  Stop;
  WriteResult(t_0,  M3x3[2], 'M3: Invert.original', SizeOf(M3x3[0]));




  for i := 0 to N-1 do    M3x3[i] := t_0;
  Start;
  for i := 0 to N-1 do
    Invert2(M3x3[i]);
  Stop;
  WriteResult(t_0,  M3x3[2], 'M3: Invert2', SizeOf(M3x3[0]));




  for i := 0 to N-1 do    M3x3[i] := t_0;
  Start;
  for i := 0 to N-1 do
    Invert_delphi_assembler_cleared(M3x3[i]);
  Stop;
  WriteResult(S1,  S1, 'M3: Invert_delphi_assembler_cleared', SizeOf(M3x3[0]));




  for i := 0 to N-1 do    M3x3[i] := t_0;
  Start;
  for i := 0 to N-1 do
    Invert_gauss(M3x3[i]);
  Stop;
  WriteResult(t_0,  M3x3[2], 'M3: Invert_gauss', SizeOf(M3x3[0]));



  for i := 0 to N-1 do    M3x3[i] := t_0;
  Start;
  for i := 0 to N-1 do
    T_SSE.Invert_gauss(M3x3[i]);
  Stop;
  WriteResult(t_0,  M3x3[2], 'M3: T_SSE.Invert_gauss', SizeOf(M3x3[0]));



  for i := 0 to N-1 do   M3x3[i] := t_0;
  Start;
    for i := 0 to N-1 do
      T_SSE.Invert(M3x3[I] );
  Stop;
  WriteResult(t_0,  M3x3[2], 'M3: T_SSE.Invert',SizeOf(M3x3[0]));




  for i := 0 to N-1 do   M3x3compact[i] := t_0;
  Start;
    for i := 0 to N-1 do
      Invert( M3x3compact[I] );
  Stop;

  WriteResult(t_0,  M3x3compact[2].ToTens, 'M3: Invert(packed matrix)', SizeOf(M3x3compact[0]));


  for i := 0 to N-1 do   M3x3[i] := t_0;
  Start;
    T_SSE.Invert(M3x3[0], 0, N );
  Stop;
  WriteResult(t_0,  M3x3[2], 'M3: T_SSE.Invert(N)', SizeOf(M3x3[0]));


  for i := 0 to N-1 do   M3x3compact[i] := t_0;
  Start;
    for i := 0 to N-1 do
      M3x3compact[I] := T_SSE.Invert_M3_clang( M3x3compact[I] );
  Stop;
  WriteResult( t_0,  M3x3compact[2].ToTens, 'M3: T_SSE.Invert_CLANG', SizeOf(M3x3compact[0]));



end;

procedure TForm1.Test_Mij_Vj;
var
    a: T_Vect;
    i: integer;
begin
  Memo1.Lines.Add('');

  zero(a);

  Start;
  for i := 0 to N-1 do    Add(a,  M3x3[i], V3[i]);
  Stop;

  WriteResult(
    Format( 'Add(a.fixed , Mij, Vj).classic: res= %9.8g; %9.8g; %9.8g;',
        [a.x, a.y, a.z]),
    SizeOf(M3x3[0]) + SizeOf(V3[0])
  );

  ////////////////////

  Start;
  zero(a);
  for i := 0 to N-1 do
    T_SSE.Add(a,  M3x3[i], V3[i]);
  Stop;

  WriteResult(
    Format( 'Add(a.fixed , Mij, Vj).SSE: res= %9.8g; %9.8g; %9.8g;',
      [a.x, a.y, a.z]),
    SizeOf(M3x3[0]) + SizeOf(V3[0])
  );

  zero(a);

  /////////////////////////////
  ///
  ///
  zero(V1);
  Start;
  for i := 0 to N-1 do    Add(V1[i],  M3x3[i], V3[i]);
  Stop;

  WriteResult( 'Add(a[k] , Mij[k], Vj[k]).classic:', SizeOf(M3x3[0]) + 2*SizeOf(V3[0]) );



  zero(V1);
  Start;
    T_SSE.Add(V1[0], M3x3[0],  V3[0],  N );
  Stop;

  WriteResult( 'T_SSE.Add(a[0] , Mij[0], Vj[0], N):', SizeOf(M3x3[0]) + 2*SizeOf(V3[0]) );

end;

procedure TForm1.Test_Vj_Mji;
var
    a: T_Vect;
    i: integer;
begin
  Memo1.Lines.Add('');

  // Test summation of V3.j*M.ji to one vector
  Start;
  zero(a);
  for i := 0 to N-1 do
    Add(a, mult(V3[i],  M3x3[i]));

  Stop;

  WriteResult(
    Format(
      'Add(a , Vj, Mji ).compiler:  res= %9.8g; %9.8g; %9.8g;',
      [a.x, a.y, a.z] ),  SizeOf(M3x3[0]) + SizeOf(V3[0])
  );


  Start;
  zero(a);
  for i := 0 to N-1 do
    T_SSE.Add(a, V3[i],  M3x3[i]);

  Stop;

  WriteResult(
    Format(
      'Add(a , Vj, Mji ).SSE:  res= %9.8g; %9.8g; %9.8g;',
      [a.x, a.y, a.z] ),  SizeOf(M3x3[0]) + SizeOf(V3[0])
  );

  zero(a);
  Start;
    T_SSE.Add(a, V3,  M3x3,  N );
  Stop;

  WriteResult(
    Format(
      'Add(a , Vj, Mji , N).SSE:  res= %9.8g; %9.8g; %9.8g;',
      [a.x, a.y, a.z] ),  SizeOf(M3x3[0]) + SizeOf(V3[0])
  );


  // Test Additon V3.j*M.ji to Array of vectors
  Start;
  zero(a);
  for i := 0 to N-1 do
    Add(V1[i], mult(V3[i],  M3x3[i]));

  Stop;

  WriteResult( 'Add(ai , Vj, Mji ).compiler:',  SizeOf(M3x3[0]) + 2*SizeOf(V3[0]) );


  Start;
  zero(a);
  for i := 0 to N-1 do
    T_SSE.Add(V1[i], V3[i],  M3x3[i]);

  Stop;

  WriteResult( 'Add(ai , Vj, Mji ).SSE:',  SizeOf(M3x3[0]) + 2*SizeOf(V3[0]) );

  zero(a);
  Start;
    T_SSE.Add(V1, V3,  M3x3,  N-1 );
  Stop;

  WriteResult( 'Add(ai , Vj, Mji, N ).SSE:',  SizeOf(M3x3[0]) + 2*SizeOf(V3[0]) );

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





end.
