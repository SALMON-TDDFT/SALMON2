# 交換計算の一時配列の寿命短縮

計算式・経路選択・MPI通信順序・継続利用するキャッシュは変更しない。
局所allocatableは手続き終了時に自動解放されるため、出口での冗長な解放追加は行わず、
最後の使用と大きな後段確保が離れている配列に限定して整理した。

- 対称性FFT再利用のforward表は実空間orbit所有者確定後、lookup表はk点置換確定後、
  inverse表はslot確定後に解放。大きな通信・密度・軌道バッファと重ねて保持しない。
- 通常の分散k点交換でもinverse表をslot作成直後に解放。
- SR wave-firstでは未使用のsource_columnを確保しない。独立したカーネルの
  コピーができた後に元multiplierを解放してからwave処理を行う。
- 通常のSR・compact経路では全格子FFT用density/spectrumを、実際の
  全格子フォールバックに初めて入る時まで確保しない。
- スクリーニングのsource_norms/global_normsはscalar envelope確定後に解放。

後段の同時保持量の削減は、default integerが4 byteの場合、対称性再利用で
`4*N_r*N_sym + 8*N_k` byte/rank、SR wave処理で`24*N_r` byte/rank。
全格子フォールバックがない通常SR/compact処理ではさらにFFT作業配列
`32*N_r*min(4,N_target)` byte/rankを確保しない。
これらは配列payloadの見積もりであり、実RSSの測定値ではない。

## 検証

記録先は `work/exchange-scratch-lifetime-20261017`。

- MPI/HSE/LibXC/ScaLAPACK、serial HSE、serial non-HSEの3ビルドが成功。
- 変更したk点交換モジュールをbounds check付きで再ビルド。
  HSE/非遮蔽カーネル、既定/明示Coulomb半径、通常/シフト付きkメッシュ、
  strided/contiguous FFT、並進操作、1/2/4 MPIの18比較が合格。
- SR directのペア再利用off/onとhaloの小規模2 MPI RTは、変更前出力との差0。
- 通常SR、全格子FFTフォールバック、occupied-supportゼロ許容スクリーニングも正常終了。
  source-supportの明示pair screeningは既存の入力仕様で拒否されたため、
  スクリーニング分岐の検証には対応するoccupied-support入力を使用。
- 422/423/424/440の計12工程が合格。
- 独立レビューで解放後参照、残存alias、未完了通信の具体的問題なし。

長時間計算・性能測定・pushは行っていない。Fujitsu・GPU・EigenExaは未検証。
