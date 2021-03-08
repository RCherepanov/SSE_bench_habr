object Form1: TForm1
  Left = 0
  Top = 0
  Caption = 'Form1'
  ClientHeight = 819
  ClientWidth = 815
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  Menu = MainMenu1
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object Memo1: TMemo
    Left = 4
    Top = 4
    Width = 809
    Height = 593
    Lines.Strings = (
      'Memo1')
    TabOrder = 0
  end
  object btn_DotProd_summ: TButton
    Left = 8
    Top = 603
    Width = 157
    Height = 25
    Caption = 'Summ of dot.prod'
    TabOrder = 1
    OnClick = btn_DotProd_summClick
  end
  object bnt_Mij_Vj: TButton
    Left = 183
    Top = 603
    Width = 157
    Height = 25
    Caption = 'Mij*Vj'
    TabOrder = 2
    OnClick = bnt_Mij_VjClick
  end
  object bnt_ViBij: TButton
    Left = 183
    Top = 628
    Width = 157
    Height = 25
    Caption = 'Vi*Mij'
    TabOrder = 3
    OnClick = bnt_ViBijClick
  end
  object btn_TestFastInvert: TButton
    Left = 0
    Top = 700
    Width = 161
    Height = 25
    Caption = 'Test fast invert 3x3'
    TabOrder = 4
    OnClick = btn_TestFastInvertClick
  end
  object btn_TestFastInvert_4x4: TButton
    Left = 183
    Top = 700
    Width = 157
    Height = 25
    Caption = 'Test fast invert 4x4'
    TabOrder = 5
    OnClick = btn_TestFastInvert_4x4Click
  end
  object btn_DotProd_Array: TButton
    Left = 4
    Top = 628
    Width = 157
    Height = 25
    Caption = 'Array Dot.prod'
    TabOrder = 6
    OnClick = btn_DotProd_ArrayClick
  end
  object MainMenu1: TMainMenu
    Left = 448
    Top = 384
    object File1: TMenuItem
      Caption = 'File'
    end
    object est1: TMenuItem
      Caption = 'Test'
    end
  end
end
