# lorenz_equation_solver
Calculation program of lorenz equation using Runge-Kutta 4 method

カオス的振る舞いが特徴のローレンツ方程式を
ルンゲクッタ４次で計算します。

定数、初期条件等はソースコード中の定数をいじることで
変更できます。

実行するとdatファイルが作成されますので、
gnuplotでプロットします。

$ gnuplot

gnuplot> splot "lorenz_eq.dat" w l
