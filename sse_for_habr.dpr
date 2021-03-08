program sse_for_habr;

uses
  FastMM4,
  Vcl.Forms,
  sse_Main in 'sse_Main.pas' {Form1},
  base.types in 'base.types.pas';

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.
