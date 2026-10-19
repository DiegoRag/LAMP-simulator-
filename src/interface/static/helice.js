/* helice.js — Modelo 3D de DNA desenhado em Canvas (B-DNA com regiões LAMP) */
(function () {
  "use strict";

  // Paletas de cor para temas claro e escuro (RGB)
  var PALETAS = {
    claro: {
      fio: [126, 144, 194],
      degrau: [126, 144, 194],
      brilho: [232, 236, 248],
      externo: [17, 47, 131],
      fip: [230, 57, 127],
      fipClaro: [242, 140, 179],
      bip: [219, 134, 14],
      bipClaro: [233, 169, 77],
      loop: [46, 125, 50],
      outro: [90, 90, 90],
      textoRotulo: {
        externo: "#112F83",
        fip: "#A3154F",
        fipClaro: "#A3154F",
        bip: "#8A5300",
        bipClaro: "#8A5300",
        loop: "#2E7D32",
        outro: "#5A5A5A"
      },
      regua: "rgba(90, 90, 90, .55)"
    },
    escuro: {
      fio: [88, 106, 168],
      degrau: [88, 106, 168],
      brilho: [196, 208, 245],
      externo: [214, 224, 255],
      fip: [255, 111, 165],
      fipClaro: [255, 160, 198],
      bip: [255, 179, 71],
      bipClaro: [255, 206, 140],
      loop: [61, 220, 151],
      outro: [180, 180, 180],
      textoRotulo: {
        externo: "#D6E0FF",
        fip: "#FF9CC2",
        fipClaro: "#FF9CC2",
        bip: "#FFC56B",
        bipClaro: "#FFC56B",
        loop: "#7BE8B8",
        outro: "#D9D9D9"
      },
      regua: "rgba(200, 210, 240, .45)"
    }
  };

  // As oito regiões do conjunto LAMP sobre a fita-alvo
  var REGIOES_MODELO = [
    { nome: "F3",  ini: 0.05, fim: 0.12, familia: "externo" },
    { nome: "F2",  ini: 0.16, fim: 0.23, familia: "fip" },
    { nome: "LF",  ini: 0.26, fim: 0.32, familia: "loop" },
    { nome: "F1",  ini: 0.35, fim: 0.42, familia: "fipClaro" },
    { nome: "B1c", ini: 0.58, fim: 0.65, familia: "bipClaro" },
    { nome: "LB",  ini: 0.68, fim: 0.74, familia: "loop" },
    { nome: "B2c", ini: 0.77, fim: 0.84, familia: "bip" },
    { nome: "B3c", ini: 0.88, fim: 0.95, familia: "externo" }
  ];

  function clamp(v, min, max) {
    return v < min ? min : (v > max ? max : v);
  }

  function suave(u) {
    var x = clamp((u - 0.55) / 0.45, 0, 1);
    return x * x * (3 - 2 * x);
  }

  function misturarRgb(c1, c2, frac) {
    var f = clamp(frac, 0, 1);
    return [
      Math.round(c1[0] + (c2[0] - c1[0]) * f),
      Math.round(c1[1] + (c2[1] - c1[1]) * f),
      Math.round(c1[2] + (c2[2] - c1[2]) * f)
    ];
  }

  // Tipos de peças para o arranjo ordenado por profundidade (desempenho A6)
  var TIPO_SEGMENTO = 1;
  var TIPO_DEGRAU = 2;
  var TIPO_ESFERA = 3;
  var TIPO_TRACEJADO = 4;

  function compararChave(a, b) {
    return a.chave - b.chave;
  }

  // Cria um controlador de animação da hélice de DNA
  function criar(canvas, opcoes) {
    if (!canvas) {
      throw new Error("Canvas é obrigatório para LampHelice.criar");
    }

    opcoes = opcoes || {};
    var tema = opcoes.tema === "escuro" ? "escuro" : "claro";
    var velocidadePadrao = tema === "escuro" ? 0.8 : 0.35;
    var velocidade = typeof opcoes.velocidade === "number" ? opcoes.velocidade : velocidadePadrao;
    var passo = typeof opcoes.passo === "number" ? opcoes.passo : 18;
    var regua = opcoes.regua !== undefined ? Boolean(opcoes.regua) : true;
    var paleta = PALETAS[tema];

    var ctx = canvas.getContext("2d");
    var W = 0;
    var H = 0;
    var phi = 0;
    var ultimoTempo = 0;
    var animacaoId = null;
    var rodando = false;
    var ativo = false;
    var emPausa = false;

    var reduzido = document.documentElement.classList.contains("movimento-reduzido");

    var comprimento = (typeof opcoes.comprimento === "number" && !isNaN(opcoes.comprimento) && opcoes.comprimento > 0)
      ? Math.round(opcoes.comprimento)
      : null;

    // Inicializa o estado de encaixe das regiões
    var listaBase;
    if (opcoes.regioes !== undefined && opcoes.regioes !== null) {
      listaBase = [];
      var fonteRegioes = Array.isArray(opcoes.regioes) ? opcoes.regioes : [];
      for (var ir = 0; ir < fonteRegioes.length; ir++) {
        var itemReg = fonteRegioes[ir];
        if (!itemReg) continue;
        var iniVal = typeof itemReg.ini === "number" ? itemReg.ini : parseFloat(itemReg.ini);
        var fimVal = typeof itemReg.fim === "number" ? itemReg.fim : parseFloat(itemReg.fim);
        if (!isNaN(iniVal) && !isNaN(fimVal) && iniVal >= 0 && fimVal <= 1 && fimVal > iniVal) {
          listaBase.push({
            nome: itemReg.nome != null ? String(itemReg.nome) : "",
            ini: iniVal,
            fim: fimVal,
            familia: itemReg.familia || "outro",
            fantasma: Boolean(itemReg.fantasma)
          });
        }
      }
    } else {
      listaBase = REGIOES_MODELO.map(function (r) {
        return {
          nome: r.nome,
          ini: r.ini,
          fim: r.fim,
          familia: r.familia,
          fantasma: false
        };
      });
    }

    var regioes = listaBase.map(function (r) {
      return {
        nome: r.nome,
        ini: r.ini,
        fim: r.fim,
        familia: r.familia,
        fantasma: Boolean(r.fantasma),
        disparado: false,
        t_inicio: 0,
        u: 0,
        ease_u: 0,
        alfa_fantasma: 0
      };
    });

    function acharRegiao(fk) {
      var candidato = null;
      for (var i = 0; i < regioes.length; i++) {
        var r = regioes[i];
        if (fk >= r.ini && fk <= r.fim) {
          if (!r.fantasma) return r;
          if (!candidato) candidato = r;
        }
      }
      return candidato;
    }

    function redimensionar() {
      var dpr = window.devicePixelRatio || 1;
      var rect = canvas.getBoundingClientRect();
      var w = rect.width || canvas.clientWidth;
      var h = rect.height || canvas.clientHeight;
      if (!w || !h) return;

      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      W = w;
      H = h;
    }

    // Arranjo de peças e buffers reaproveitados entre quadros (desempenho A6)
    var pecas = [];
    var numPecas = 0;
    var amostrasA = [];
    var amostrasB = [];
    var bufferPrimer = [];

    function obterPeca() {
      if (numPecas < pecas.length) {
        var p = pecas[numPecas];
        p.pontos = null;
        numPecas++;
        return p;
      }
      var nova = {
        tipo: 0,
        chave: 0,
        x1: 0,
        y1: 0,
        x2: 0,
        y2: 0,
        x: 0,
        y: 0,
        r: 0,
        lw: 0,
        alfa: 0,
        cr: 0,
        cg: 0,
        cb: 0,
        pontos: null
      };
      pecas.push(nova);
      numPecas++;
      return nova;
    }

    // Desenha um quadro individual da hélice
    function desenharQuadro(agora) {
      if (!agora) agora = performance.now();
      if (W <= 0 || H <= 0) return;

      ctx.clearRect(0, 0, W, H);

      var mx = Math.max(24, W * 0.03);
      var x0 = mx;
      var x1 = W - mx;
      var cy = H * 0.5;
      var R = Math.min(H * 0.24, 78);

      var N = clamp(Math.round((x1 - x0) / passo), 24, 140);
      var F = 900;

      // Linha do eixo, tracejada, em cy
      ctx.save();
      ctx.strokeStyle = "rgba(" + paleta.fio.join(",") + ", 0.08)";
      ctx.lineWidth = 1;
      ctx.setLineDash([2, 6]);
      ctx.beginPath();
      ctx.moveTo(x0, cy);
      ctx.lineTo(x1, cy);
      ctx.stroke();
      ctx.restore();

      // Régua em pb
      if (regua) {
        var yRegua = cy + R + 30;

        ctx.save();
        ctx.strokeStyle = paleta.regua;
        ctx.fillStyle = paleta.regua;
        ctx.lineWidth = 1;
        ctx.font = "700 10.5px 'Courier New', monospace";
        ctx.textAlign = "center";
        ctx.textBaseline = "top";

        if (comprimento) {
          var PASSOS = [10, 20, 25, 40, 50, 100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000];
          var passoPb = PASSOS[PASSOS.length - 1];
          for (var pi = 0; pi < PASSOS.length; pi++) {
            if (comprimento / PASSOS[pi] <= 8) {
              passoPb = PASSOS[pi];
              break;
            }
          }

          var valoresPb = [];
          var vPb = 0;
          while (vPb < comprimento) {
            valoresPb.push(vPb);
            vPb += passoPb;
          }
          valoresPb.push(comprimento);

          var distMarcas = valoresPb.length > 1
            ? (x1 - x0) * (valoresPb[1] - valoresPb[0]) / comprimento
            : (x1 - x0);

          var passoPular = 1;
          if (distMarcas < 46) {
            passoPular = 2;
            if (distMarcas * 2 < 46) {
              passoPular = 4;
            }
          }

          var ultimoIdx = valoresPb.length - 1;
          var penultimoIdx = -1;
          for (var j = ultimoIdx - 1; j >= 0; j--) {
            if (j % passoPular === 0) {
              penultimoIdx = j;
              break;
            }
          }

          var pularPenultimo = false;
          if (penultimoIdx >= 0) {
            var xUltimo = x0 + (x1 - x0) * valoresPb[ultimoIdx] / comprimento;
            var xPenultimo = x0 + (x1 - x0) * valoresPb[penultimoIdx] / comprimento;
            var distPenultimo = xUltimo - xPenultimo;

            var textoUltimo = valoresPb[ultimoIdx] + " pb";
            var textoPenultimo = String(valoresPb[penultimoIdx]);
            var wUltimo = ctx.measureText(textoUltimo).width;
            var wPenultimo = ctx.measureText(textoPenultimo).width;

            if (distPenultimo < 46 || distPenultimo < (wPenultimo + wUltimo) / 2) {
              pularPenultimo = true;
            }
          }

          for (var ki = 0; ki < valoresPb.length; ki++) {
            var valPb = valoresPb[ki];
            var xkPb = x0 + (x1 - x0) * valPb / comprimento;
            ctx.beginPath();
            ctx.moveTo(xkPb, yRegua);
            ctx.lineTo(xkPb, yRegua + 7);
            ctx.stroke();

            var ehUltimo = (ki === ultimoIdx);
            var ehPenultimoPulado = (ki === penultimoIdx && pularPenultimo);
            var deveNumerar = ehUltimo || (!ehPenultimoPulado && (ki % passoPular === 0));

            if (deveNumerar) {
              // coordenadas a partir de 1, como no GenBank: o início do alvo é a base 1
              var textoPb = ehUltimo ? valPb + " pb" : String(valPb === 0 ? 1 : valPb);
              ctx.measureText(textoPb);
              ctx.fillText(textoPb, xkPb, yRegua + 10);
            }
          }
        } else {
          var dxR = (x1 - x0) * 20 / (N - 1);
          var passoPularR = 1;
          if (dxR < 46) {
            passoPularR = 2;
            if (dxR * 2 < 46) {
              passoPularR = 4;
            }
          }
          var intervaloR = 20 * passoPularR;

          var ultR = -1;
          for (var idx = 0; idx < N; idx++) {
            if (idx % intervaloR === 0) ultR = idx;
          }

          var penultR = ultR >= intervaloR ? ultR - intervaloR : -1;
          var pularPenultR = false;
          if (penultR >= 0) {
            var xkUltR = x0 + (x1 - x0) * ultR / (N - 1);
            var xkPenultR = x0 + (x1 - x0) * penultR / (N - 1);
            var distPenultR = xkUltR - xkPenultR;
            var textoUltR = ultR + " pb";
            var textoPenultR = String(penultR);
            var wUltR = ctx.measureText(textoUltR).width;
            var wPenultR = ctx.measureText(textoPenultR).width;
            if (distPenultR < 46 || distPenultR < (wPenultR + wUltR) / 2) {
              pularPenultR = true;
            }
          }

          for (var kR = 0; kR < N; kR++) {
            if (kR % 10 === 0) {
              var xkR = x0 + (x1 - x0) * kR / (N - 1);
              var ehUltR = (kR === ultR);
              var ehPenultPuladoR = (kR === penultR && pularPenultR);
              var deveNumerarR = ehUltR || (!ehPenultPuladoR && (kR % intervaloR === 0));

              if (kR % 20 === 0) {
                ctx.beginPath();
                ctx.moveTo(xkR, yRegua);
                ctx.lineTo(xkR, yRegua + 7);
                ctx.stroke();
                if (deveNumerarR) {
                  var textoNum = ehUltR ? kR + " pb" : String(kR);
                  ctx.measureText(textoNum);
                  ctx.fillText(textoNum, xkR, yRegua + 10);
                }
              } else {
                ctx.beginPath();
                ctx.moveTo(xkR, yRegua);
                ctx.lineTo(xkR, yRegua + 4);
                ctx.stroke();
                if (deveNumerarR) {
                  var textoNum2 = ehUltR ? kR + " pb" : String(kR);
                  ctx.measureText(textoNum2);
                  ctx.fillText(textoNum2, xkR, yRegua + 10);
                }
              }
            }
          }
        }
        ctx.restore();
      }

      // Atualiza o progresso de encaixe de cada região
      for (var rIdx = 0; rIdx < regioes.length; rIdx++) {
        var reg = regioes[rIdx];
        if (reg.disparado) {
          var decorrido = agora - reg.t_inicio;
          if (reg.fantasma) {
            reg.u = 0;
            reg.ease_u = 0;
            if (decorrido < 0) {
              reg.alfa_fantasma = 0;
            } else if (decorrido >= 600) {
              reg.alfa_fantasma = 0.55;
            } else {
              reg.alfa_fantasma = 0.55 * (decorrido / 600);
            }
          } else {
            if (decorrido < 0) {
              reg.u = 0;
            } else if (decorrido >= 900) {
              reg.u = 1;
            } else {
              reg.u = decorrido / 900;
            }
            reg.ease_u = 1 - Math.pow(1 - reg.u, 3);
          }
        } else {
          reg.u = 0;
          reg.ease_u = 0;
          reg.alfa_fantasma = 0;
        }
      }

      // A1. Fitas lisas: amostragem com M = 6 subpassos por par
      var M = 6;
      var totalAmostras = (N - 1) * M + 1;

      for (var j = 0; j < totalAmostras; j++) {
        var u = j / M;
        var x = x0 + (x1 - x0) * u / (N - 1);
        var f = u / (N - 1);
        var e = clamp(Math.min((x - x0) / 110, (x1 - x) / 110), 0, 1);
        var theta = 2 * Math.PI * u / 10.5 + phi;

        // Fita A em θ
        var Y_A = R * Math.sin(theta);
        var Z_A = R * Math.cos(theta);
        var s_A = F / (F - Z_A);
        var sx_A = x;
        var sy_A = cy + Y_A * s_A;
        var t_A = (Z_A / R + 1) / 2;

        if (!amostrasA[j]) amostrasA[j] = {};
        var aA = amostrasA[j];
        aA.j = j;
        aA.u = u;
        aA.x = x;
        aA.f = f;
        aA.e = e;
        aA.theta = theta;
        aA.Y = Y_A;
        aA.Z = Z_A;
        aA.s = s_A;
        aA.sx = sx_A;
        aA.sy = sy_A;
        aA.t = t_A;

        // Fita B em θ + 0.78π
        var theta_B = theta + 0.78 * Math.PI;
        var Y_B = R * Math.sin(theta_B);
        var Z_B = R * Math.cos(theta_B);
        var s_B = F / (F - Z_B);
        var sx_B = x;
        var sy_B = cy + Y_B * s_B;
        var t_B = (Z_B / R + 1) / 2;

        if (!amostrasB[j]) amostrasB[j] = {};
        var aB = amostrasB[j];
        aB.j = j;
        aB.u = u;
        aB.x = x;
        aB.f = f;
        aB.e = e;
        aB.theta = theta_B;
        aB.Y = Y_B;
        aB.Z = Z_B;
        aB.s = s_B;
        aB.sx = sx_B;
        aB.sy = sy_B;
        aB.t = t_B;
      }

      numPecas = 0;

      // Segmentos da fita A
      for (var jA = 0; jA < totalAmostras - 1; jA++) {
        var p1A = amostrasA[jA];
        var p2A = amostrasA[jA + 1];
        var t_mA = (p1A.t + p2A.t) / 2;
        var s_mA = (p1A.s + p2A.s) / 2;
        var e_mA = (p1A.e + p2A.e) / 2;
        var f_mA = (p1A.f + p2A.f) / 2;

        var regA = acharRegiao(f_mA);
        var corA = paleta.fio;
        if (regA && !regA.fantasma && regA.ease_u > 0) {
          corA = misturarRgb(paleta.fio, paleta[regA.familia] || paleta.fio, regA.ease_u);
        }

        var pecaA = obterPeca();
        pecaA.tipo = TIPO_SEGMENTO;
        pecaA.chave = (p1A.Z + p2A.Z) / 2;
        pecaA.x1 = p1A.sx;
        pecaA.y1 = p1A.sy;
        pecaA.x2 = p2A.sx;
        pecaA.y2 = p2A.sy;
        pecaA.lw = (1.6 + 2.8 * t_mA) * s_mA;
        pecaA.alfa = (0.18 + 0.72 * t_mA) * e_mA;
        pecaA.cr = corA[0];
        pecaA.cg = corA[1];
        pecaA.cb = corA[2];
      }

      // Segmentos da fita B (sempre cor do fio)
      for (var jB = 0; jB < totalAmostras - 1; jB++) {
        var p1B = amostrasB[jB];
        var p2B = amostrasB[jB + 1];
        var t_mB = (p1B.t + p2B.t) / 2;
        var s_mB = (p1B.s + p2B.s) / 2;
        var e_mB = (p1B.e + p2B.e) / 2;
        var corB = paleta.fio;

        var pecaB = obterPeca();
        pecaB.tipo = TIPO_SEGMENTO;
        pecaB.chave = (p1B.Z + p2B.Z) / 2;
        pecaB.x1 = p1B.sx;
        pecaB.y1 = p1B.sy;
        pecaB.x2 = p2B.sx;
        pecaB.y2 = p2B.sy;
        pecaB.lw = (1.6 + 2.8 * t_mB) * s_mB;
        pecaB.alfa = (0.18 + 0.72 * t_mB) * e_mB;
        pecaB.cr = corB[0];
        pecaB.cg = corB[1];
        pecaB.cb = corB[2];
      }

      // A2. Degraus em duas metades (pares inteiros k = 0 … N - 1)
      for (var kD = 0; kD < N; kD++) {
        var pAk = amostrasA[kD * M];
        var pBk = amostrasB[kD * M];

        var Z_MD = (pAk.Z + pBk.Z) / 2;
        var Y_MD = (pAk.Y + pBk.Y) / 2;
        var s_MD = F / (F - Z_MD);
        var sx_MD = pAk.x;
        var sy_MD = cy + Y_MD * s_MD;

        // Metade A: de A_k até M_k com a luz da fita A
        var t_DA = pAk.t;
        var s_DA = pAk.s;
        var e_DA = pAk.e;
        var regDA = acharRegiao(pAk.f);
        var corDegA = paleta.degrau;
        var alfaDegA = (0.08 + 0.34 * t_DA) * e_DA;

        if (regDA && !regDA.fantasma && regDA.ease_u > 0) {
          var euDA = regDA.ease_u;
          corDegA = misturarRgb(paleta.degrau, paleta[regDA.familia] || paleta.degrau, euDA);
          var alfaNormA = (0.08 + 0.34 * t_DA) * e_DA;
          var alfaFamA = (0.2 + 0.5 * t_DA) * e_DA;
          alfaDegA = alfaNormA + (alfaFamA - alfaNormA) * euDA;
        }

        var pecaDA = obterPeca();
        pecaDA.tipo = TIPO_DEGRAU;
        pecaDA.chave = (pAk.Z + Z_MD) / 2 - 0.5;
        pecaDA.x1 = pAk.sx;
        pecaDA.y1 = pAk.sy;
        pecaDA.x2 = sx_MD;
        pecaDA.y2 = sy_MD;
        pecaDA.lw = (1.4 + 1.4 * t_DA) * s_DA;
        pecaDA.alfa = alfaDegA;
        pecaDA.cr = corDegA[0];
        pecaDA.cg = corDegA[1];
        pecaDA.cb = corDegA[2];

        // Metade B: de M_k até B_k com a luz da fita B
        var t_DB = pBk.t;
        var s_DB = pBk.s;
        var e_DB = pBk.e;
        var corDegB = paleta.degrau;
        var alfaDegB = (0.08 + 0.34 * t_DB) * e_DB;

        var pecaDB = obterPeca();
        pecaDB.tipo = TIPO_DEGRAU;
        pecaDB.chave = (pBk.Z + Z_MD) / 2 - 0.5;
        pecaDB.x1 = sx_MD;
        pecaDB.y1 = sy_MD;
        pecaDB.x2 = pBk.sx;
        pecaDB.y2 = pBk.sy;
        pecaDB.lw = (1.4 + 1.4 * t_DB) * s_DB;
        pecaDB.alfa = alfaDegB;
        pecaDB.cr = corDegB[0];
        pecaDB.cg = corDegB[1];
        pecaDB.cb = corDegB[2];
      }

      // A3. Nucleotídeos (esferas só nos pares inteiros k = 0 … N - 1)
      for (var kN = 0; kN < N; kN++) {
        // Fita A
        var ptNA = amostrasA[kN * M];
        var rNA = (2.2 + 2.6 * ptNA.t) * ptNA.s;
        var alfaNA = (0.30 + 0.70 * ptNA.t) * ptNA.e;
        var regNA = acharRegiao(ptNA.f);
        var corNA = paleta.fio;
        if (regNA && !regNA.fantasma && regNA.ease_u > 0) {
          corNA = misturarRgb(paleta.fio, paleta[regNA.familia] || paleta.fio, regNA.ease_u);
        }

        var pecaNA = obterPeca();
        pecaNA.tipo = TIPO_ESFERA;
        pecaNA.chave = ptNA.Z;
        pecaNA.x = ptNA.sx;
        pecaNA.y = ptNA.sy;
        pecaNA.r = rNA;
        pecaNA.alfa = alfaNA;
        pecaNA.cr = corNA[0];
        pecaNA.cg = corNA[1];
        pecaNA.cb = corNA[2];

        // Fita B
        var ptNB = amostrasB[kN * M];
        var rNB = (2.2 + 2.6 * ptNB.t) * ptNB.s;
        var alfaNB = (0.30 + 0.70 * ptNB.t) * ptNB.e;
        var corNB = paleta.fio;

        var pecaNB = obterPeca();
        pecaNB.tipo = TIPO_ESFERA;
        pecaNB.chave = ptNB.Z;
        pecaNB.x = ptNB.sx;
        pecaNB.y = ptNB.sy;
        pecaNB.r = rNB;
        pecaNB.alfa = alfaNB;
        pecaNB.cr = corNB[0];
        pecaNB.cg = corNB[1];
        pecaNB.cb = corNB[2];
      }

      // A4. Primer que encaixa (hibridiza radialmente enquanto 0 < u < 1) ou flutua fantasma
      for (var rP = 0; rP < regioes.length; rP++) {
        var regPrimer = regioes[rP];
        var desenharNormal = !regPrimer.fantasma && regPrimer.u > 0 && regPrimer.u < 1;
        var desenharFantasma = regPrimer.fantasma && regPrimer.alfa_fantasma > 0;

        if (desenharNormal || desenharFantasma) {
          var easeP = regPrimer.fantasma ? 0 : regPrimer.ease_u;
          var ru = R * (1 + 0.32 * (1 - easeP));
          var alfaPrimer = regPrimer.fantasma ? regPrimer.alfa_fantasma : Math.min(1, 2 * regPrimer.u);
          var corPrimer = paleta[regPrimer.familia] || paleta.fio;

          var numPontosPrimer = 0;
          for (var jP = 0; jP < totalAmostras; jP++) {
            var ptP = amostrasA[jP];
            if (ptP.f >= regPrimer.ini - 1e-5 && ptP.f <= regPrimer.fim + 1e-5) {
              var Y_P = ru * Math.sin(ptP.theta);
              var Z_P = ru * Math.cos(ptP.theta);
              var s_P = F / (F - Z_P);
              var sx_P = ptP.x;
              var sy_P = cy + Y_P * s_P;
              var t_P = (Z_P / ru + 1) / 2;

              if (numPontosPrimer < bufferPrimer.length) {
                var itemP = bufferPrimer[numPontosPrimer];
                itemP.px = sx_P;
                itemP.py = sy_P;
                itemP.t = t_P;
                itemP.s = s_P;
                itemP.Z = Z_P;
                itemP.isInteger = (ptP.j % M === 0);
              } else {
                bufferPrimer.push({
                  px: sx_P,
                  py: sy_P,
                  t: t_P,
                  s: s_P,
                  Z: Z_P,
                  isInteger: (ptP.j % M === 0)
                });
              }
              numPontosPrimer++;
            }
          }

          if (regPrimer.fantasma) {
            // Linha tracejada única contínua com setLineDash([5, 5])
            if (numPontosPrimer > 1) {
              var zTotal = 0;
              var sTotal = 0;
              var tTotal = 0;
              var pontosTrac = [];
              for (var pi = 0; pi < numPontosPrimer; pi++) {
                var pTrac = bufferPrimer[pi];
                zTotal += pTrac.Z;
                sTotal += pTrac.s;
                tTotal += pTrac.t;
                pontosTrac.push({ x: pTrac.px, y: pTrac.py });
              }
              var tMedio = tTotal / numPontosPrimer;
              var sMedio = sTotal / numPontosPrimer;
              var chaveTrac = (zTotal / numPontosPrimer) + 300 * (1 - easeP);

              var pecaTrac = obterPeca();
              pecaTrac.tipo = TIPO_TRACEJADO;
              pecaTrac.chave = chaveTrac;
              pecaTrac.pontos = pontosTrac;
              pecaTrac.lw = (2.4 + 2.4 * tMedio) * sMedio;
              pecaTrac.alfa = alfaPrimer;
              pecaTrac.cr = corPrimer[0];
              pecaTrac.cg = corPrimer[1];
              pecaTrac.cb = corPrimer[2];
            }

            // Esferas nos pares inteiros com metade do alfa
            for (var ePf = 0; ePf < numPontosPrimer; ePf++) {
              var pntPf = bufferPrimer[ePf];
              if (pntPf.isInteger) {
                var chaveEsfPf = pntPf.Z + 300 * (1 - easeP);
                var rEsfPf = (2.4 + 2.6 * pntPf.t) * pntPf.s;

                var pecaEsfPf = obterPeca();
                pecaEsfPf.tipo = TIPO_ESFERA;
                pecaEsfPf.chave = chaveEsfPf;
                pecaEsfPf.x = pntPf.px;
                pecaEsfPf.y = pntPf.py;
                pecaEsfPf.r = rEsfPf;
                pecaEsfPf.alfa = alfaPrimer * 0.5;
                pecaEsfPf.cr = corPrimer[0];
                pecaEsfPf.cg = corPrimer[1];
                pecaEsfPf.cb = corPrimer[2];
              }
            }
          } else {
            // Linha lisa por essas posições
            for (var sP = 0; sP < numPontosPrimer - 1; sP++) {
              var p1P = bufferPrimer[sP];
              var p2P = bufferPrimer[sP + 1];
              var t_mP = (p1P.t + p2P.t) / 2;
              var s_mP = (p1P.s + p2P.s) / 2;
              var chaveSegP = (p1P.Z + p2P.Z) / 2 + 300 * (1 - easeP);

              var pecaSegP = obterPeca();
              pecaSegP.tipo = TIPO_SEGMENTO;
              pecaSegP.chave = chaveSegP;
              pecaSegP.x1 = p1P.px;
              pecaSegP.y1 = p1P.py;
              pecaSegP.x2 = p2P.px;
              pecaSegP.y2 = p2P.py;
              pecaSegP.lw = (2.4 + 2.4 * t_mP) * s_mP;
              pecaSegP.alfa = alfaPrimer;
              pecaSegP.cr = corPrimer[0];
              pecaSegP.cg = corPrimer[1];
              pecaSegP.cb = corPrimer[2];
            }

            // Esferas nos pares inteiros
            for (var eP = 0; eP < numPontosPrimer; eP++) {
              var pntP = bufferPrimer[eP];
              if (pntP.isInteger) {
                var chaveEsfP = pntP.Z + 300 * (1 - easeP);
                var rEsfP = (2.4 + 2.6 * pntP.t) * pntP.s;

                var pecaEsfP = obterPeca();
                pecaEsfP.tipo = TIPO_ESFERA;
                pecaEsfP.chave = chaveEsfP;
                pecaEsfP.x = pntP.px;
                pecaEsfP.y = pntP.py;
                pecaEsfP.r = rEsfP;
                pecaEsfP.alfa = alfaPrimer;
                pecaEsfP.cr = corPrimer[0];
                pecaEsfP.cg = corPrimer[1];
                pecaEsfP.cb = corPrimer[2];
              }
            }
          }
        }
      }

      // Ordena e desenha as peças pelo pintor (mais fundo primeiro)
      pecas.length = numPecas;
      pecas.sort(compararChave);

      for (var pIdx = 0; pIdx < numPecas; pIdx++) {
        var peca = pecas[pIdx];
        switch (peca.tipo) {
          case TIPO_SEGMENTO:
            ctx.lineWidth = peca.lw;
            ctx.strokeStyle = "rgba(" + peca.cr + "," + peca.cg + "," + peca.cb + "," + peca.alfa + ")";
            ctx.lineCap = "round";
            ctx.beginPath();
            ctx.moveTo(peca.x1, peca.y1);
            ctx.lineTo(peca.x2, peca.y2);
            ctx.stroke();
            break;

          case TIPO_DEGRAU:
            ctx.lineWidth = peca.lw;
            ctx.strokeStyle = "rgba(" + peca.cr + "," + peca.cg + "," + peca.cb + "," + peca.alfa + ")";
            ctx.lineCap = "butt";
            ctx.beginPath();
            ctx.moveTo(peca.x1, peca.y1);
            ctx.lineTo(peca.x2, peca.y2);
            ctx.stroke();
            break;

          case TIPO_ESFERA:
            if (peca.alfa > 0 && peca.r > 0) {
              var gx = peca.x - 0.35 * peca.r;
              var gy = peca.y - 0.35 * peca.r;
              var grad = ctx.createRadialGradient(gx, gy, 0.1 * peca.r, peca.x, peca.y, peca.r);
              grad.addColorStop(0, "rgba(" + paleta.brilho.join(",") + "," + peca.alfa + ")");
              grad.addColorStop(0.55, "rgba(" + peca.cr + "," + peca.cg + "," + peca.cb + "," + peca.alfa + ")");
              grad.addColorStop(1, "rgba(" + Math.round(peca.cr * 0.7) + "," +
                Math.round(peca.cg * 0.7) + "," + Math.round(peca.cb * 0.7) + "," + peca.alfa + ")");
              ctx.fillStyle = grad;
              ctx.beginPath();
              ctx.arc(peca.x, peca.y, peca.r, 0, 2 * Math.PI);
              ctx.fill();
            }
            break;

          case TIPO_TRACEJADO:
            if (peca.alfa > 0 && peca.pontos && peca.pontos.length > 1) {
              ctx.save();
              ctx.lineWidth = peca.lw;
              ctx.strokeStyle = "rgba(" + peca.cr + "," + peca.cg + "," + peca.cb + "," + peca.alfa + ")";
              ctx.lineCap = "round";
              ctx.setLineDash([5, 5]);
              ctx.beginPath();
              ctx.moveTo(peca.pontos[0].x, peca.pontos[0].y);
              for (var ptIdx = 1; ptIdx < peca.pontos.length; ptIdx++) {
                ctx.lineTo(peca.pontos[ptIdx].x, peca.pontos[ptIdx].y);
              }
              ctx.stroke();
              ctx.setLineDash([]);
              ctx.restore();
            }
            break;
        }
      }

      // Anel de encaixe ao completar (600 ms após 900 ms)
      for (var rH = 0; rH < regioes.length; rH++) {
        var regAnel = regioes[rH];
        if (regAnel.disparado && !regAnel.fantasma) {
          var tempoDecorrido = agora - regAnel.t_inicio;
          if (tempoDecorrido >= 900 && tempoDecorrido <= 1500) {
            var vAnel = (tempoDecorrido - 900) / 600;
            var alfaAnel = 0.5 * (1 - vAnel);
            var fCentroH = (regAnel.ini + regAnel.fim) / 2;
            var sx_cH = x0 + (x1 - x0) * fCentroH;
            var sy_cH = cy;
            var wH = (x1 - x0) * (regAnel.fim - regAnel.ini);
            var rxH = wH * (0.55 + 0.40 * vAnel);

            var uCentroH = fCentroH * (N - 1);
            var thetaCentroH = 2 * Math.PI * uCentroH / 10.5 + phi;
            var Z_cH = R * Math.cos(thetaCentroH);
            var s_cH = F / (F - Z_cH);
            var ryH = R * s_cH * (0.95 + 0.50 * vAnel);

            var lwH = 2 - 1 * vAnel;
            var corAnel = paleta[regAnel.familia] || paleta.fio;

            ctx.save();
            ctx.strokeStyle = "rgba(" + corAnel.join(",") + "," + alfaAnel + ")";
            ctx.lineWidth = lwH;
            ctx.beginPath();
            ctx.ellipse(sx_cH, sy_cH, rxH, ryH, 0, 0, 2 * Math.PI);
            ctx.stroke();
            ctx.restore();
          }
        }
      }

      // Rótulos, colchetes e fios das regiões encaixadas
      var deltaLinhas = (cy - R - 46 - 8 < 0) ? 12 : 16;
      var sy_topoL = cy - R;
      var sy_colcheteL = cy - (R + 16);
      var sy_rotulo_baixo = cy - (R + 30);
      var sy_rotulo_cima = cy - (R + 30 + deltaLinhas);

      var itensRotulo = [];
      for (var rL = 0; rL < regioes.length; rL++) {
        var regRot = regioes[rL];
        var alfaRot;
        if (regRot.fantasma) {
          if (regRot.disparado) {
            var decRot = agora - regRot.t_inicio;
            if (decRot >= 600) {
              alfaRot = 1;
            } else if (decRot >= 0) {
              alfaRot = decRot / 600;
            } else {
              alfaRot = 0;
            }
          } else {
            alfaRot = 0;
          }
        } else {
          alfaRot = suave(regRot.u);
        }

        if (alfaRot > 0) {
          var fCentroL = (regRot.ini + regRot.fim) / 2;
          var sx_cL = x0 + (x1 - x0) * fCentroL;
          var sx_iniL = x0 + (x1 - x0) * regRot.ini;
          var sx_fimL = x0 + (x1 - x0) * regRot.fim;

          var textoNome = regRot.nome;
          var comecaSem = /^sem /i.test(textoNome);
          if (regRot.fantasma && !comecaSem && textoNome.slice(-1) !== "?") {
            textoNome = textoNome + "?";
          }

          var fonteRotulo = regRot.fantasma
            ? "italic 700 12.5px 'Open Sans', Arial, sans-serif"
            : (tema === "escuro"
                ? "700 13px 'Open Sans', Arial, sans-serif"
                : "700 12.5px 'Open Sans', Arial, sans-serif");

          ctx.font = fonteRotulo;
          var largRotulo = ctx.measureText(textoNome).width;

          itensRotulo.push({
            reg: regRot,
            alfa: alfaRot,
            indice: rL,
            sx_c: sx_cL,
            sx_ini: sx_iniL,
            sx_fim: sx_fimL,
            texto: textoNome,
            fonte: fonteRotulo,
            largura: largRotulo,
            ladoEsq: sx_cL - largRotulo / 2,
            ladoDir: sx_cL + largRotulo / 2,
            linha: "baixo"
          });
        }
      }

      itensRotulo.sort(function (a, b) {
        return a.sx_c !== b.sx_c ? a.sx_c - b.sx_c : a.indice - b.indice;
      });

      var ultimoDirBaixo = null;
      var ultimoDirCima = null;

      for (var iR = 0; iR < itensRotulo.length; iR++) {
        var itemRot = itensRotulo[iR];
        var cabeBaixo = (ultimoDirBaixo === null) || (itemRot.ladoEsq >= ultimoDirBaixo + 6);
        if (cabeBaixo) {
          itemRot.linha = "baixo";
          ultimoDirBaixo = itemRot.ladoDir;
        } else {
          var cabeCima = (ultimoDirCima === null) || (itemRot.ladoEsq >= ultimoDirCima + 6);
          if (cabeCima) {
            itemRot.linha = "cima";
            ultimoDirCima = itemRot.ladoDir;
          } else {
            itemRot.linha = "baixo";
            ultimoDirBaixo = Math.max(ultimoDirBaixo, itemRot.ladoDir);
          }
        }
      }

      for (var iD = 0; iD < itensRotulo.length; iD++) {
        var itemDes = itensRotulo[iD];
        var ehCima = itemDes.linha === "cima";
        var sy_rot = ehCima ? sy_rotulo_cima : sy_rotulo_baixo;
        var sy_fio = ehCima ? (sy_colcheteL - deltaLinhas) : sy_colcheteL;
        var corRot = paleta[itemDes.reg.familia] || paleta.fio;

        // Fio de 1px com alfa 0.4
        ctx.save();
        ctx.strokeStyle = "rgba(" + corRot.join(",") + "," + (0.4 * itemDes.alfa) + ")";
        ctx.lineWidth = 1;
        if (itemDes.reg.fantasma) {
          ctx.setLineDash([5, 5]);
        }
        ctx.beginPath();
        ctx.moveTo(itemDes.sx_c, sy_fio);
        ctx.lineTo(itemDes.sx_c, sy_topoL);
        ctx.stroke();
        if (itemDes.reg.fantasma) {
          ctx.setLineDash([]);
        }
        ctx.restore();

        // Colchete de 2px
        ctx.save();
        ctx.strokeStyle = "rgba(" + corRot.join(",") + "," + itemDes.alfa + ")";
        ctx.lineWidth = 2;
        ctx.lineCap = "round";
        if (itemDes.reg.fantasma) {
          ctx.setLineDash([5, 5]);
        }
        ctx.beginPath();
        ctx.moveTo(itemDes.sx_ini, sy_colcheteL + 4);
        ctx.lineTo(itemDes.sx_ini, sy_colcheteL);
        ctx.lineTo(itemDes.sx_fim, sy_colcheteL);
        ctx.lineTo(itemDes.sx_fim, sy_colcheteL + 4);
        ctx.stroke();
        if (itemDes.reg.fantasma) {
          ctx.setLineDash([]);
        }
        ctx.restore();

        // Rótulo de texto centralizado acima da hélice
        var corTexto = (paleta.textoRotulo && paleta.textoRotulo[itemDes.reg.familia])
          ? paleta.textoRotulo[itemDes.reg.familia]
          : "rgb(" + corRot.join(",") + ")";

        ctx.save();
        ctx.globalAlpha = itemDes.alfa;
        ctx.font = itemDes.fonte;
        ctx.fillStyle = corTexto;
        ctx.textAlign = "center";
        ctx.textBaseline = "middle";
        ctx.fillText(itemDes.texto, itemDes.sx_c, sy_rot);
        ctx.restore();
      }
    }

    // Laço principal com requestAnimationFrame
    function laco(tempoAtual) {
      if (!rodando) return;

      var dt = ultimoTempo ? (tempoAtual - ultimoTempo) / 1000 : 0.016;
      ultimoTempo = tempoAtual;
      if (dt > 0.05) dt = 0.05;

      phi += velocidade * dt;
      desenharQuadro(tempoAtual);

      animacaoId = requestAnimationFrame(laco);
    }

    function iniciar() {
      if (reduzido) {
        phi = 0.6;
        for (var i = 0; i < regioes.length; i++) {
          regioes[i].disparado = true;
          if (!regioes[i].fantasma) {
            regioes[i].u = 1;
            regioes[i].ease_u = 1;
          } else {
            regioes[i].u = 0;
            regioes[i].ease_u = 0;
            regioes[i].alfa_fantasma = 0.55;
          }
          regioes[i].t_inicio = performance.now() - 2000;
        }
        redimensionar();
        desenharQuadro(performance.now());
        return;
      }

      ativo = true;
      if (rodando) return;

      redimensionar();
      rodando = true;
      ultimoTempo = performance.now();
      animacaoId = requestAnimationFrame(laco);
    }

    function parar() {
      rodando = false;
      ativo = false;
      if (animacaoId) {
        cancelAnimationFrame(animacaoId);
        animacaoId = null;
      }
    }

    function encaixarEmSequencia(intervaloMs) {
      var intervalo = typeof intervaloMs === "number" ? intervaloMs : 170;
      var agora = performance.now();
      for (var i = 0; i < regioes.length; i++) {
        regioes[i].disparado = true;
        regioes[i].t_inicio = agora + i * intervalo;
      }
      if (!rodando) {
        desenharQuadro(agora);
      }
    }

    function definirProgresso(p, imediato) {
      var progresso = clamp(p, 0, 1);
      var limite = Math.floor(progresso * regioes.length + 1e-6);
      var agora = performance.now();
      for (var i = 0; i < regioes.length; i++) {
        if (regioes[i].fantasma) {
          continue;
        }
        if (i < limite) {
          if (imediato === true) {
            regioes[i].disparado = true;
            regioes[i].u = 1;
            regioes[i].ease_u = 1;
            regioes[i].t_inicio = agora - 2000;
          } else if (!regioes[i].disparado) {
            regioes[i].disparado = true;
            regioes[i].t_inicio = agora;
          }
        } else if (imediato === true) {
          regioes[i].disparado = false;
          regioes[i].u = 0;
          regioes[i].ease_u = 0;
          regioes[i].t_inicio = 0;
        }
      }
      if (!rodando) {
        desenharQuadro(agora);
      }
    }

    function encaixarTudo() {
      var agora = performance.now();
      for (var i = 0; i < regioes.length; i++) {
        regioes[i].disparado = true;
        if (!regioes[i].fantasma) {
          regioes[i].u = 1;
          regioes[i].ease_u = 1;
        } else {
          regioes[i].u = 0;
          regioes[i].ease_u = 0;
          regioes[i].alfa_fantasma = 0.55;
        }
        regioes[i].t_inicio = agora - 2000;
      }
      if (!rodando) {
        desenharQuadro(agora);
      }
    }

    // Observadores de tamanho e visibilidade
    var ro = null;
    if (typeof ResizeObserver !== "undefined") {
      ro = new ResizeObserver(function () {
        redimensionar();
        if (!rodando) {
          desenharQuadro(performance.now());
        }
      });
      ro.observe(canvas);
    }

    var io = null;
    if (typeof IntersectionObserver !== "undefined") {
      io = new IntersectionObserver(function (entries) {
        entries.forEach(function (entry) {
          if (entry.isIntersecting) {
            if (ativo && emPausa) {
              emPausa = false;
              rodando = true;
              ultimoTempo = performance.now();
              animacaoId = requestAnimationFrame(laco);
            }
          } else {
            if (rodando) {
              emPausa = true;
              rodando = false;
              if (animacaoId) {
                cancelAnimationFrame(animacaoId);
                animacaoId = null;
              }
            }
          }
        });
      }, { threshold: 0.01 });
      io.observe(canvas);
    }

    function onVisibilityChange() {
      if (document.hidden) {
        if (rodando) {
          emPausa = true;
          rodando = false;
          if (animacaoId) {
            cancelAnimationFrame(animacaoId);
            animacaoId = null;
          }
        }
      } else {
        if (ativo && emPausa) {
          emPausa = false;
          rodando = true;
          ultimoTempo = performance.now();
          animacaoId = requestAnimationFrame(laco);
        }
      }
    }

    document.addEventListener("visibilitychange", onVisibilityChange);

    function destruir() {
      parar();
      if (ro) {
        ro.disconnect();
        ro = null;
      }
      if (io) {
        io.disconnect();
        io = null;
      }
      document.removeEventListener("visibilitychange", onVisibilityChange);
    }

    // Inicialização da geometria inicial
    redimensionar();
    if (reduzido) {
      phi = 0.6;
      encaixarTudo();
    } else {
      desenharQuadro(performance.now());
    }

    return {
      iniciar: iniciar,
      parar: parar,
      encaixarEmSequencia: encaixarEmSequencia,
      definirProgresso: definirProgresso,
      encaixarTudo: encaixarTudo,
      destruir: destruir
    };
  }

  window.LampHelice = {
    criar: criar
  };
})();
