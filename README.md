
# NOisy Read Assembler (nora) とは

オックスフォードナノポア用のデノボアセンブラです。<br>
アキュラシーが70～80% 程度の長いリードをインプットとして想定しており、コンティグをFASTA形式で出力します。<br>
オーバーラップ・レイアウト・コンセンサス法に分類されるアセンブラです。<br>

# E. coli のデータで環状コンティグになった

![contigs](./img/ont_pruned.png)
相補鎖を消す前のコンティグの様子 ([Bandage](https://github.com/rrwick/Bandage) を使用)

![dotplot1](./img/blastdotplot.png)
相補鎖を消した後のコンティグと参照ゲノムとのドットプロット([blast2dotplot.py](https://qiita.com/satoshi_kawato/items/0e5c13621e53bad8d9a0) を使用)

![dotplot2](./img/2,c_U00096.resized.png)
相補鎖を消した後のコンティグと参照ゲノムとのドットプロット([moddotplot](https://github.com/marbl/ModDotPlot) を使用)

[Cali et al](https://pubmed.ncbi.nlm.nih.gov/29617724/) にナノポアのリードが載っていました、
[Loman labs](http://lab.loman.net/2016/07/30/nanopore-r9-data-release/) の下の方のR9 1D Rapid FASTA files (1.5Gb)というデータです。<br>
nora でアセンブルしたコンティグをMG1655とドットプロットしたら上のようになりました（ドットプロットの見方は [A quick reference guide for interpretting the dot plot](https://mummer.sourceforge.net/manual/#mummerplot)の' here'のリンク先のpdfが分かりやすいと思います）。<br>
コンティグ一本のサイズは4.4 Mb なのでMG1655の4.6 Mb より短いですが、Cali et al のTable 10 を見るとMiniasm もこのくらいなので、まずまずといったところでしょうか。<br>


# 今後
今後はE. coli 以外のデータもアセンブルしてみようと思います。<br>
そして使いやすいアセンブラを作りつつ、系統樹など？応用にも挑戦してみたいです。<br>

# 方法

今回のE. coli をアセンブルした方法は<br>
`time ./bin/nora.pl --prefix foo ont_data.fa --genome_size 5000000 --bin_dir ./bin/`<br>
です（ont\_data.fa は Loman labs からダウンロードしてきたFASTA ファイルです。また ./bin/ にバイナリを配置してあります）。<br>
foo.contigs.fa に２本のコンティグが出たので、コンティグ同士のドットプロットを参考に相補鎖と判断し１本を選びました。<br>

実行環境は<br>
CPU: AMD Ryzen 5 3500 6-Core Processor 3.59 GHz (6 core 6 threads)<br>
RAM: 64 GB<br>
Windows 11, WSL2 上のUbuntu<br>
で76 min (Wall clock time), 397 min (CPU time) かかりました。

## インストール方法

`git clone https://github.com/inco-imai/nora.git`<br>
`cd nora`<br>
`mkdir bin`<br>
`make`<br>
`make install`<br>
CPUはリトルエンディアンを前提としています。<br>
RAMは64GB以上にして下さい（32GBでも動くかもしれませんが試していません）<br>

# フロー
![nora's flow](./img/flow.png)

詳しくはnora.pl を読んで下さい。

# 各ファイル解説

## mhseol (Min Hash Seed and Extend OverLapper: mhseol\_v3.c, b\_heikki.dynamic.h, mymurmurhash3.h, myqsort.h)
mhseolは[MHAP](https://www.nature.com/articles/nbt.3238) に倣ってMinHashを使っています。MinHashは簡単に言うと二本のリード間で共有されるkmerの個数が高速に分かるフィルターみたいなものです。<br>
ハッシュ関数は[murmurhash3](https://en.wikipedia.org/wiki/MurmurHash)で散らした後[XORShift RNG](https://en.wikipedia.org/wiki/Xorshift)を使っています（この高速化のアイディアもMHAP同様です）。<br>

Seed and Extend は[BLAST](https://pubmed.ncbi.nlm.nih.gov/2231712/)に倣っています。<br>
Extendでアラインメントする所は[Heikki 2002](http://www.stringology.org/event/2002/p6.html)を実装しました。このアルゴリズムは簡単に言うとビットパラレルなバンディッドアラインメントです。
また、CHUNK\_SIZE(デフォルトは100)ごとに一致率を見て閾値未満なら相同性検索を打ち切り、打ち切らない場合はバンドの中心をずらす実装にしてあります（これは自分で考えました）。<br>
トレースバックももちろんします。これも自分で考えました。<br>

バンドの中心をずらしつつバンディッドアラインメントするイメージ図
![dynamic banded heikki](./img/dynamic_banded_heikki.png)

途中で並列クイックソートをしますが、これは[stanaka2](https://qiita.com/stanaka2/items/526102c6c56759c22b8f)さんに倣いました。OpenMPを使った並列クイックソートです。

mhseolは下図のようにクエリqがサブジェクトsに含まれてしまう場合のアラインメントは出力しません。<br>
また-Oオプションを外せばqがsを完全に含む場合（上から３つめ）も出力しないよう選べます。<br>
![mhseol\_pass\_and\_discard](./img/mhseol_pass_and_discard.png)


## error correction (bfmt72s, nss2v, vertical2cfq\_and\_blacklist.pl)
![error correction](./img/ec.png)
noraのエラーコレクションは、図のようにmhseolで作ったペアワイズアラインメントを集めてきて縦に積み多数決を取ることで行います。<br>

ブラックリストblにはキメラが疑われるリードやデプスがmin\_depth未満の領域があるリードを入れます。グラフトラバース時に使われません。<br>
![blacklist](./img/blacklist.png)


## gt (Graph Traverser: gt\_v3.c)
noraはリードをノードとし、オーバーハングのあるオーバーラップをエッジとしてグラフを作ります。<br>
またエラーコレクションの時に作ったblに入っているノードを無効にします。<br>

その後バブルポッピング(bp)を下図のように行います。<br>
![bubble popping](./img/bp.png)
noraはデフォルトではまずD=2,3,4を行い、グラフを整理してからD=2,3,4,...,6\*4=24 まで行います。こうすることでいきなりD=24まで行うより速くなりました。<br>
ちなみに24という数字に根拠はありません。<br>

また、下図のようなループは削除します。bp同様にD=6\*4=24まで行いました(こちらではD=1も行います）。<br>
![loop removal](./img/loop_removal.png)

bp後は枝刈りpruningを行います、つまりエッジが分岐している場合、その分岐しているノードをスタートとしてエッジ数が最大のルートを残し、他のルートはlimit(30)以下のエッジ数なら消します。
それでもまだ分岐が残っている場合、分岐のあるノードをRAM上のbl（≠ファイルのbl）に加えて、グラフを元に戻した後再度bp→pruningを行います。
blが変わらなくなる（収束する）か10回繰り返すとbp→pruningは終わります。

### 視覚化
gt\_v3 と後の処理でprefix.orig.cd.gfa, prefix.init.cd.gfa, prefix.bp.cd.gfa, prefix.pruned.cd.gfa という[Bandage](https://github.com/rrwick/Bandage) で絵にできる４つのファイルが出力されます。
それぞれについて、origはオーバーラップグラフの初期状態を、initはブラックリストのノードを除いた状態を、bpはバブルポッピング後の状態を、prunedは枝刈りが終わった後の状態を表します。

これらの視覚化はアセンブリがうまくいかない時、どの状態に問題があるかの予想に役立ちます。
例えばorigでグラフがこま切れの時はオーバーラップをパラメータを変更してやり直す必要があると思います。
initで問題があるならmin\_depthを変えてブラックリストを作り直す必要があるかもしれません。


## その他

### bfmt7viewer.pl
bfmt7(blast format 7)のアラインメントの様子を見るコードです。lessと組み合わせて使います。

### cfq2fa.pl
cfq (consensus fastq)からエラーコレクションしたリードを取り出します。おまけ。

### clean\_gfa.pl
gfa形式からエッジのないノードを除きます。cd.gfa=cleaned gfa

### extract\_gfa.pl
相補鎖を消したgfaを作るときに使います。おまけ。

### fa2gfaS.pl
FASTA形式のリードをノードとしてgfaの一部を作ります。後にgfalとcatで組み合わせます。

### fa2idfa\_v2.pl
FASTAのリード名を数字のIDにします(idfa形式)。通常相補鎖も出力します。-fオプションをつけると相補鎖を出力しません。

### get\_blastout\_in\_xml.pl
blast2dotplot.py (外部スクリプト）のインプットを作るショートカットです。

### get\_fasta\_stats.pl
FASTAファイルのリード数、塩基数、N50などを見ます。

### get\_longest\_20x\_fa.pl
-gオプションでゲノムサイズを与え、入力のFASTAファイルで長い方のリードから20x分集めて出力します。

### idfa2samheader.pl
nss形式(name sorted sam形式)に付けるヘッダをidfa形式から作ります。おまけ。

### nora.pl
noraアセンブラのインターフェイスです。

### ohf2fa\_v2.pl
グラフトラバースしてオーバーハングする領域を集めたohf(over hanging format)フォーマットとエラーコレクションしたリードcfq(corrected fastq)を合わせてコンティグを作ります。

### vertical2fa.pl
.vertical をfastaにします。おまけ。

### vertical.vim
vimで.vertical を見るときに色を付けることでちょっとしたゲノムブラウザにできます。vim使いは ~/.vim/syntax/ フォルダにコピーしましょう。おまけ。


# contact
email: `inco.imai0@gmail.com` <br>
日本語か英語でお願いします。
