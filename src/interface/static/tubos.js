"use strict";

// Cor dos tubos ao longo do tempo (página /tubos), em minutos de reação.
// Modo "simulacao": curvas do lampsim para o conjunto analisado (só leitura).
// Modo "manual": curvas digitadas à mão, com a cor do mesmo modelo do lampsim (paleta).
(function () {
  function iniciar() {
    var bandejaEl = document.querySelector("[data-tubos-bandeja]");
    var scriptEl = document.getElementById("tubos-dados");
    if (!bandejaEl || !scriptEl) return;

    var dados;
    try {
      dados = JSON.parse(scriptEl.textContent || "{}");
    } catch (e) {
      return;
    }
    if (!dados || !Array.isArray(dados.tubos) || dados.tubos.length === 0) return;

    var MODO = dados.modo === "simulacao" ? "simulacao" : "manual";
    var CHAVE_STORAGE = "lamp-tubos-v2";
    var PASSOS_EIXO = [1, 2, 5, 10, 15, 20, 30, 60, 120];
    var PH_MIN = 6.0;
    var PH_MAX = 9.0;
    var DURACAO_MAX = numeroPositivo(dados.duracao_max_min, 60);
    var PASSO = numeroPositivo(dados.passo_min, 0.5);
    var LEITURA = numeroPositivo(dados.leitura_min, 15);
    var PKA = numeroPositivo(dados.pka, 8);
    var T_MIN = Array.isArray(dados.t_min) ? dados.t_min : [];
    var PALETA = Array.isArray(dados.paleta) && dados.paleta.length > 1 ? dados.paleta : ["#FE00B1", "#F3DF48"];

    function numeroPositivo(v, padrao) {
      return typeof v === "number" && isFinite(v) && v > 0 ? v : padrao;
    }

    function virgula(v, casas) {
      return Number(v).toFixed(casas).replace(".", ",");
    }

    // -------------------------------------------------------------
    // Cores
    // -------------------------------------------------------------

    function hexParaRgb(hex) {
      var n = parseInt(String(hex || "").replace("#", ""), 16);
      if (isNaN(n)) return [0, 0, 0];
      return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
    }

    function rgbParaHex(rgb) {
      return "#" + ((1 << 24) + (rgb[0] << 16) + (rgb[1] << 8) + rgb[2]).toString(16).slice(1).toUpperCase();
    }

    var cacheRgb = {};
    function rgbDe(hex) {
      if (!cacheRgb[hex]) cacheRgb[hex] = hexParaRgb(hex);
      return cacheRgb[hex];
    }

    function clonarManuais(lista) {
      return lista.map(function (tb) {
        return {
          nome: String(tb.nome || "Tubo"),
          curva: tb.curva === "sigmoide" ? "sigmoide" : "linear",
          inicio: Math.max(0, Number(tb.inicio) || 0),
          duracao: Number(tb.duracao) > 0 ? Number(tb.duracao) : 1,
          alcance: typeof tb.alcance === "number" ? Math.max(0, Math.min(100, tb.alcance)) : 100
        };
      });
    }

    // Progresso de uma curva manual (0 a 1) no minuto t
    function progressoManual(tubo, t) {
      var x = Math.max(0, Math.min(1, (t - tubo.inicio) / (tubo.duracao > 0 ? tubo.duracao : 0.001)));
      var f = x;
      if (tubo.curva === "sigmoide") {
        var l0 = 1 / (1 + Math.exp(5));
        var l1 = 1 / (1 + Math.exp(-5));
        f = (1 / (1 + Math.exp(-10 * (x - 0.5))) - l0) / (l1 - l0);
      }
      return Math.max(0, Math.min(1, (f * tubo.alcance) / 100));
    }

    // Valor interpolado de uma série amostrada em T_MIN (passo uniforme)
    function interpolar(serie, t) {
      if (!serie || serie.length === 0) return 0;
      var pos = Math.max(0, Math.min(serie.length - 1, t / PASSO));
      var i = Math.floor(pos);
      var j = Math.min(serie.length - 1, i + 1);
      var w = pos - i;
      return serie[i] * (1 - w) + serie[j] * w;
    }

    function corRgb(i, t) {
      var tb = estado.tubos[i];
      if (MODO === "simulacao") {
        var cores = tb.cores || [];
        if (cores.length === 0) return [0, 0, 0];
        var pos = Math.max(0, Math.min(cores.length - 1, t / PASSO));
        var a = Math.floor(pos);
        var b = Math.min(cores.length - 1, a + 1);
        var w = pos - a;
        var ca = rgbDe(cores[a]);
        var cb = rgbDe(cores[b]);
        return [
          Math.round(ca[0] * (1 - w) + cb[0] * w),
          Math.round(ca[1] * (1 - w) + cb[1] * w),
          Math.round(ca[2] * (1 - w) + cb[2] * w)
        ];
      }
      var p = progressoManual(tb, t);
      return rgbDe(PALETA[Math.round(p * (PALETA.length - 1))]);
    }

    // Valor mostrado no gráfico: pH (simulação) ou progresso 0–1 (manual)
    function valor(i, t) {
      var tb = estado.tubos[i];
      return MODO === "simulacao" ? interpolar(tb.ph, t) : progressoManual(tb, t);
    }

    function rotuloTubo(i, t) {
      if (MODO === "simulacao") return "pH " + virgula(valor(i, t), 1);
      return Math.round(valor(i, t) * 100) + " %";
    }

    // -------------------------------------------------------------
    // Estado e armazenamento local
    // -------------------------------------------------------------

    var estado = {
      duracao: Math.min(numeroPositivo(dados.duracao_min, 30), DURACAO_MAX),
      velocidade: numeroPositivo(dados.velocidade, 60),
      tubos: MODO === "simulacao" ? dados.tubos : clonarManuais(dados.tubos),
      selecionado: 0,
      t: 0
    };

    function lerSalvo() {
      try {
        var bruto = localStorage.getItem(CHAVE_STORAGE);
        var obj = bruto ? JSON.parse(bruto) : null;
        return obj && typeof obj === "object" ? obj : {};
      } catch (e) {
        return {};
      }
    }

    (function restaurar() {
      var salvo = lerSalvo()[MODO];
      if (!salvo || typeof salvo !== "object") return;
      if (typeof salvo.duracao === "number" && salvo.duracao >= 1 && salvo.duracao <= DURACAO_MAX) {
        estado.duracao = salvo.duracao;
      }
      if (typeof salvo.velocidade === "number" && salvo.velocidade > 0) estado.velocidade = salvo.velocidade;
      if (typeof salvo.selecionado === "number" && salvo.selecionado >= 0 && salvo.selecionado < estado.tubos.length) {
        estado.selecionado = Math.floor(salvo.selecionado);
      }
      if (typeof salvo.t === "number" && salvo.t >= 0) estado.t = Math.min(salvo.t, estado.duracao);
      if (MODO === "manual" && Array.isArray(salvo.tubos) && salvo.tubos.length === estado.tubos.length) {
        estado.tubos = clonarManuais(salvo.tubos);
      }
    })();

    var timerSalvar = null;
    function agendarSalvar() {
      if (timerSalvar) clearTimeout(timerSalvar);
      timerSalvar = setTimeout(function () {
        timerSalvar = null;
        var tudo = lerSalvo();
        tudo[MODO] = {
          duracao: estado.duracao,
          velocidade: estado.velocidade,
          selecionado: estado.selecionado,
          t: estado.t,
          tubos: MODO === "manual" ? estado.tubos : undefined
        };
        try {
          localStorage.setItem(CHAVE_STORAGE, JSON.stringify(tudo));
        } catch (e) {
          // sem armazenamento local
        }
      }, 200);
    }

    // -------------------------------------------------------------
    // Elementos
    // -------------------------------------------------------------

    var tubosBotoes = Array.prototype.slice.call(bandejaEl.querySelectorAll(".tubo"));
    var tocarBotao = document.querySelector("[data-tubos-tocar]");
    var cursorInput = document.querySelector("[data-tubos-cursor]");
    var relogioOutput = document.querySelector("[data-tubos-relogio]");
    var duracaoInput = document.querySelector("[data-tubos-duracao]");
    var velocidadeSelect = document.querySelector("[data-tubos-velocidade]");
    var faixasLinhas = Array.prototype.slice.call(document.querySelectorAll(".tubos-faixa-linha"));
    var faixasCanvases = Array.prototype.slice.call(document.querySelectorAll("[data-tubos-faixa]"));
    var marcadorEl = document.querySelector("[data-tubos-marcador]");
    var marcaLeituraEl = document.querySelector("[data-tubos-marca-leitura]");
    var eixoEl = document.querySelector("[data-tubos-eixo]");
    var linhasTabela = Array.prototype.slice.call(document.querySelectorAll("[data-tubos-linha]"));
    var editorEl = document.querySelector("[data-tubos-editor]");
    var editorNomeEl = document.querySelector("[data-tubos-editor-nome]");
    var curvaCanvas = document.querySelector("[data-tubos-curva]");
    var leituraEl = document.querySelector("[data-tubos-leitura]");
    var copiarBotao = document.querySelector("[data-tubos-copiar]");
    var restaurarBotao = document.querySelector("[data-tubos-restaurar]");

    var faixasSujas = estado.tubos.map(function () { return true; });
    function sujarFaixas() {
      for (var i = 0; i < faixasSujas.length; i++) faixasSujas[i] = true;
    }

    // -------------------------------------------------------------
    // Eixo e marcadores das faixas
    // -------------------------------------------------------------

    function recriarEixo() {
      if (!eixoEl) return;
      eixoEl.innerHTML = "";
      var largura = eixoEl.clientWidth || (faixasCanvases[0] ? faixasCanvases[0].clientWidth : 0);
      var alvo = Math.max(3, Math.min(7, Math.floor(largura / 64)));
      var passo = PASSOS_EIXO[PASSOS_EIXO.length - 1];
      for (var k = 0; k < PASSOS_EIXO.length; k++) {
        if (Math.floor(estado.duracao / PASSOS_EIXO[k]) + 1 <= alvo) {
          passo = PASSOS_EIXO[k];
          break;
        }
      }
      for (var m = 0; m <= estado.duracao + 1e-9; m += passo) {
        var span = document.createElement("span");
        span.style.setProperty("--x", (m / estado.duracao).toFixed(5));
        span.textContent = m + " min";
        eixoEl.appendChild(span);
      }
      if (marcaLeituraEl) {
        var dentro = LEITURA <= estado.duracao;
        marcaLeituraEl.hidden = !dentro;
        marcaLeituraEl.style.setProperty("--x", (Math.min(LEITURA, estado.duracao) / estado.duracao).toFixed(5));
      }
    }

    // -------------------------------------------------------------
    // Desenho
    // -------------------------------------------------------------

    var offscreen = document.createElement("canvas");
    var offCtx = offscreen.getContext("2d");

    function faixaDeCores(i, largura) {
      if (offscreen.width < largura) offscreen.width = largura;
      offscreen.height = 1;
      var img = offCtx.createImageData(largura, 1);
      for (var col = 0; col < largura; col++) {
        var rgb = corRgb(i, (col / Math.max(1, largura - 1)) * estado.duracao);
        img.data[col * 4] = rgb[0];
        img.data[col * 4 + 1] = rgb[1];
        img.data[col * 4 + 2] = rgb[2];
        img.data[col * 4 + 3] = 255;
      }
      offCtx.putImageData(img, 0, 0);
    }

    function renderizarFaixas() {
      var dpr = window.devicePixelRatio || 1;
      for (var i = 0; i < faixasCanvases.length; i++) {
        if (!faixasSujas[i] || !estado.tubos[i]) continue;
        var canvas = faixasCanvases[i];
        if (canvas.clientWidth <= 0 || canvas.clientHeight <= 0) continue;
        var w = Math.max(1, Math.round(canvas.clientWidth * dpr));
        var h = Math.max(1, Math.round(canvas.clientHeight * dpr));
        if (canvas.width !== w || canvas.height !== h) {
          canvas.width = w;
          canvas.height = h;
        }
        faixaDeCores(i, w);
        var ctx = canvas.getContext("2d");
        ctx.imageSmoothingEnabled = false;
        ctx.drawImage(offscreen, 0, 0, w, 1, 0, 0, w, h);
        faixasSujas[i] = false;
      }
    }

    function renderizarCurva() {
      if (!curvaCanvas || curvaCanvas.clientWidth <= 0 || curvaCanvas.clientHeight <= 0) return;
      var dpr = window.devicePixelRatio || 1;
      var wCss = curvaCanvas.clientWidth;
      var hCss = curvaCanvas.clientHeight;
      if (curvaCanvas.width !== Math.round(wCss * dpr) || curvaCanvas.height !== Math.round(hCss * dpr)) {
        curvaCanvas.width = Math.round(wCss * dpr);
        curvaCanvas.height = Math.round(hCss * dpr);
      }
      var ctx = curvaCanvas.getContext("2d");
      ctx.save();
      ctx.scale(dpr, dpr);
      ctx.clearRect(0, 0, wCss, hCss);

      var x0 = 52;
      var x1 = Math.max(x0 + 10, wCss - 16);
      var y0 = 14;
      var y1 = Math.max(y0 + 10, hCss - 34);
      var pw = x1 - x0;
      var ph = y1 - y0;
      var i = estado.selecionado;
      var dur = estado.duracao;

      function xDe(t) { return x0 + (t / dur) * pw; }
      function yDe(v) {
        if (MODO === "simulacao") return y0 + ((PH_MAX - Math.max(PH_MIN, Math.min(PH_MAX, v))) / (PH_MAX - PH_MIN)) * ph;
        return y1 - v * ph;
      }

      // grade e rótulos do eixo vertical
      ctx.font = '12px "Open Sans", sans-serif';
      ctx.fillStyle = "#5A5A5A";
      ctx.strokeStyle = "#D9D9D9";
      ctx.lineWidth = 1;
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      var marcas = MODO === "simulacao"
        ? [[9, "pH 9"], [PKA, "pKa " + virgula(PKA, 1)], [7, "pH 7"], [6, "pH 6"]]
        : [[1, "100 %"], [0.5, "50 %"], [0, "0 %"]];
      ctx.save();
      ctx.setLineDash([4, 4]);
      marcas.forEach(function (m) {
        var y = yDe(m[0]);
        ctx.beginPath();
        ctx.moveTo(x0, y);
        ctx.lineTo(x1, y);
        ctx.stroke();
        ctx.fillText(m[1], x0 - 8, y);
      });
      ctx.restore();

      // leitura nominal (simulação)
      if (MODO === "simulacao" && LEITURA <= dur) {
        ctx.save();
        ctx.strokeStyle = "#112F83";
        ctx.setLineDash([6, 4]);
        ctx.beginPath();
        ctx.moveTo(xDe(LEITURA), y0);
        ctx.lineTo(xDe(LEITURA), y1);
        ctx.stroke();
        ctx.restore();
      }

      // eixo do tempo
      ctx.textBaseline = "top";
      ctx.textAlign = "left";
      ctx.fillText("0 min", x0, y1 + 18);
      ctx.textAlign = "right";
      ctx.fillText(virgula(dur, dur % 1 === 0 ? 0 : 1) + " min", x1, y1 + 18);

      // faixa de cor sob o gráfico
      var passos = Math.max(10, Math.ceil(pw));
      faixaDeCores(i, passos);
      ctx.imageSmoothingEnabled = false;
      ctx.drawImage(offscreen, 0, 0, passos, 1, x0, y1 + 4, pw, 10);

      // curva
      ctx.save();
      ctx.strokeStyle = "#112F83";
      ctx.lineWidth = 3;
      ctx.lineJoin = "round";
      ctx.lineCap = "round";
      ctx.beginPath();
      for (var s = 0; s <= passos; s++) {
        var t = (s / passos) * dur;
        if (s === 0) ctx.moveTo(xDe(t), yDe(valor(i, t)));
        else ctx.lineTo(xDe(t), yDe(valor(i, t)));
      }
      ctx.stroke();
      ctx.restore();

      // cursor
      var tc = Math.max(0, Math.min(dur, estado.t));
      ctx.save();
      ctx.strokeStyle = "#DB860E";
      ctx.fillStyle = "#DB860E";
      ctx.lineWidth = 1.8;
      ctx.beginPath();
      ctx.moveTo(xDe(tc), y0);
      ctx.lineTo(xDe(tc), y1);
      ctx.stroke();
      ctx.beginPath();
      ctx.arc(xDe(tc), yDe(valor(i, tc)), 5, 0, Math.PI * 2);
      ctx.fill();
      ctx.restore();

      ctx.restore();
    }

    var ultimaLeitura = 0;
    var mensagem = null;
    var timerMensagem = null;

    function atualizarLeitura(forcado) {
      if (!leituraEl) return;
      var agora = performance.now();
      if (tocando && !forcado && agora - ultimaLeitura < 250) return;
      ultimaLeitura = agora;
      if (mensagem) {
        leituraEl.textContent = mensagem;
        return;
      }
      var tb = estado.tubos[estado.selecionado];
      if (!tb) return;
      var tStr = virgula(estado.t, 1) + " min";
      if (MODO === "simulacao") {
        var vira = typeof tb.t_positivo_min === "number"
          ? "Vira em " + virgula(tb.t_positivo_min, 1) + " min"
          : "Não vira até " + virgula(DURACAO_MAX, 0) + " min";
        leituraEl.textContent = "Em " + tStr + ", " + tb.nome + " está em pH " +
          virgula(valor(estado.selecionado, estado.t), 2) + ". " + vira + "; diagnóstico aos " +
          virgula(LEITURA, 0) + " min: " + (tb.diagnostico_texto || tb.diagnostico) + ".";
        return;
      }
      if (tb.alcance === 0) {
        leituraEl.textContent = tb.nome + " não muda de cor.";
        return;
      }
      var fim = tb.inicio + tb.duracao;
      leituraEl.textContent = "Em " + tStr + ", " + tb.nome + " está a " +
        Math.round(progressoManual(tb, estado.t) * 100) + " % da queda de pH. Começa a virar em " +
        virgula(tb.inicio, 1) + " min e " +
        (fim > estado.duracao ? "termina depois do fim da janela" : "termina em " + virgula(fim, 1) + " min") + ".";
    }

    var agendado = false;
    function solicitarRender() {
      if (agendado) return;
      agendado = true;
      requestAnimationFrame(function () {
        agendado = false;
        renderizar();
      });
    }

    function renderizar() {
      for (var i = 0; i < tubosBotoes.length; i++) {
        var btn = tubosBotoes[i];
        var tb = estado.tubos[i];
        if (!tb) continue;
        var liquido = btn.querySelector(".tubo-liquido");
        if (liquido) liquido.style.setProperty("--cor", rgbParaHex(corRgb(i, estado.t)));
        var pctEl = btn.querySelector("[data-tubo-pct]");
        if (pctEl) pctEl.textContent = rotuloTubo(i, estado.t);
        var nomeEl = btn.querySelector("[data-tubo-nome]");
        if (nomeEl && nomeEl.textContent !== tb.nome) nomeEl.textContent = tb.nome;
        btn.setAttribute("aria-pressed", i === estado.selecionado ? "true" : "false");
      }
      faixasLinhas.forEach(function (linha, f) {
        linha.classList.toggle("selecionada", f === estado.selecionado);
        var tbF = estado.tubos[f];
        var nomeF = linha.querySelector("[data-faixa-nome]");
        if (tbF && nomeF && nomeF.textContent !== tbF.nome) {
          nomeF.textContent = tbF.nome;
          linha.setAttribute("aria-label", "Escolher " + tbF.nome);
        }
      });
      linhasTabela.forEach(function (linha, l) {
        linha.classList.toggle("selecionada", l === estado.selecionado);
      });
      if (marcadorEl) marcadorEl.style.setProperty("--x", (estado.t / estado.duracao).toFixed(5));
      if (relogioOutput) relogioOutput.textContent = virgula(estado.t, 1) + " min";
      if (cursorInput && document.activeElement !== cursorInput) cursorInput.value = String(estado.t);
      renderizarFaixas();
      renderizarCurva();
      atualizarLeitura(false);
    }

    // -------------------------------------------------------------
    // Seleção e editor (manual)
    // -------------------------------------------------------------

    function atualizarEditor() {
      if (!editorEl) return;
      var tb = estado.tubos[estado.selecionado];
      if (!tb) return;
      if (editorNomeEl) editorNomeEl.textContent = tb.nome;
      var nome = editorEl.querySelector('input[data-campo="nome"]');
      if (nome && document.activeElement !== nome) nome.value = tb.nome;
      editorEl.querySelectorAll('input[data-campo="curva"]').forEach(function (r) {
        r.checked = r.value === tb.curva;
      });
      ["inicio", "duracao", "alcance"].forEach(function (campo) {
        editorEl.querySelectorAll('[data-campo="' + campo + '"]').forEach(function (inp) {
          if (inp.type === "range" && campo !== "alcance") inp.max = String(estado.duracao);
          if (document.activeElement !== inp) inp.value = tb[campo];
        });
      });
    }

    function selecionar(i) {
      if (i < 0 || i >= estado.tubos.length) return;
      estado.selecionado = i;
      atualizarEditor();
      solicitarRender();
      atualizarLeitura(true);
      agendarSalvar();
    }

    tubosBotoes.forEach(function (btn, i) {
      btn.addEventListener("click", function () { selecionar(i); });
    });
    linhasTabela.forEach(function (linha, i) {
      linha.addEventListener("click", function () { selecionar(i); });
    });
    bandejaEl.addEventListener("keydown", function (e) {
      if (e.key !== "ArrowLeft" && e.key !== "ArrowRight") return;
      e.preventDefault();
      var n = estado.tubos.length;
      var novo = (estado.selecionado + (e.key === "ArrowLeft" ? -1 : 1) + n) % n;
      selecionar(novo);
      if (tubosBotoes[novo]) tubosBotoes[novo].focus();
    });

    if (editorEl) {
      var aplicar = function (alvo) {
        var campo = alvo.getAttribute("data-campo");
        var tb = estado.tubos[estado.selecionado];
        if (!campo || !tb) return;
        if (campo === "nome") {
          tb.nome = alvo.value;
          if (editorNomeEl) editorNomeEl.textContent = tb.nome || "Tubo " + (estado.selecionado + 1);
        } else if (campo === "curva") {
          tb.curva = alvo.value === "sigmoide" ? "sigmoide" : "linear";
        } else {
          var v = parseFloat(alvo.value);
          if (isNaN(v)) return;
          if (campo === "inicio") v = Math.max(0, v);
          if (campo === "duracao") v = Math.max(0.5, v);
          if (campo === "alcance") v = Math.max(0, Math.min(100, Math.round(v)));
          tb[campo] = v;
          editorEl.querySelectorAll('[data-campo="' + campo + '"]').forEach(function (par) {
            if (par !== alvo) par.value = v;
          });
        }
        faixasSujas[estado.selecionado] = true;
        solicitarRender();
        atualizarLeitura(true);
        agendarSalvar();
      };
      editorEl.addEventListener("input", function (e) { aplicar(e.target); });
      editorEl.addEventListener("change", function (e) { aplicar(e.target); });
      var nomeInput = editorEl.querySelector('input[data-campo="nome"]');
      if (nomeInput) {
        nomeInput.addEventListener("blur", function () {
          var tb = estado.tubos[estado.selecionado];
          if (tb && !nomeInput.value.trim()) {
            tb.nome = "Tubo " + (estado.selecionado + 1);
            nomeInput.value = tb.nome;
            if (editorNomeEl) editorNomeEl.textContent = tb.nome;
            solicitarRender();
            agendarSalvar();
          }
        });
      }
    }

    if (copiarBotao) {
      copiarBotao.addEventListener("click", function () {
        var sel = estado.tubos[estado.selecionado];
        if (!sel) return;
        estado.tubos.forEach(function (tb) {
          tb.curva = sel.curva;
          tb.inicio = sel.inicio;
          tb.duracao = sel.duracao;
          tb.alcance = sel.alcance;
        });
        sujarFaixas();
        mensagem = "Copiado para os " + estado.tubos.length + " tubos.";
        if (timerMensagem) clearTimeout(timerMensagem);
        timerMensagem = setTimeout(function () {
          mensagem = null;
          atualizarLeitura(true);
        }, 3000);
        solicitarRender();
        atualizarLeitura(true);
        agendarSalvar();
      });
    }

    if (restaurarBotao) {
      restaurarBotao.addEventListener("click", function () {
        parar();
        var tudo = lerSalvo();
        delete tudo[MODO];
        try {
          localStorage.setItem(CHAVE_STORAGE, JSON.stringify(tudo));
        } catch (e) {
          // sem armazenamento local
        }
        estado.tubos = clonarManuais(dados.tubos);
        estado.duracao = Math.min(numeroPositivo(dados.duracao_min, 30), DURACAO_MAX);
        estado.selecionado = 0;
        estado.t = 0;
        aplicarDuracao();
        atualizarEditor();
        atualizarLeitura(true);
      });
    }

    // -------------------------------------------------------------
    // Tempo: cursor, duração, velocidade e reprodução
    // -------------------------------------------------------------

    var tocando = false;
    var rafId = null;
    var anterior = 0;

    function parar() {
      if (!tocando) return;
      tocando = false;
      if (rafId !== null) cancelAnimationFrame(rafId);
      rafId = null;
      if (tocarBotao) {
        tocarBotao.setAttribute("aria-pressed", "false");
        tocarBotao.textContent = "Reproduzir";
      }
      solicitarRender();
      atualizarLeitura(true);
      agendarSalvar();
    }

    function tocar() {
      if (tocando) return;
      if (estado.t >= estado.duracao) estado.t = 0;
      tocando = true;
      if (tocarBotao) {
        tocarBotao.setAttribute("aria-pressed", "true");
        tocarBotao.textContent = "Pausar";
      }
      anterior = performance.now();
      function passo(agora) {
        if (!tocando) return;
        // velocidade: segundos de reação por segundo de tela; estado.t está em minutos
        estado.t += ((agora - anterior) / 1000) * (estado.velocidade / 60);
        anterior = agora;
        if (estado.t >= estado.duracao) {
          estado.t = estado.duracao;
          parar();
          return;
        }
        solicitarRender();
        rafId = requestAnimationFrame(passo);
      }
      rafId = requestAnimationFrame(passo);
    }

    if (tocarBotao) {
      tocarBotao.addEventListener("click", function () {
        if (tocando) parar();
        else tocar();
      });
    }

    function aplicarDuracao() {
      if (estado.t > estado.duracao) estado.t = estado.duracao;
      if (duracaoInput) duracaoInput.value = String(estado.duracao);
      if (cursorInput) {
        cursorInput.max = String(estado.duracao);
        cursorInput.step = String(PASSO);
        cursorInput.value = String(estado.t);
      }
      recriarEixo();
      sujarFaixas();
      atualizarEditor();
      solicitarRender();
    }

    if (cursorInput) {
      cursorInput.addEventListener("input", function () {
        if (tocando) parar();
        var v = parseFloat(cursorInput.value);
        if (isNaN(v)) return;
        estado.t = Math.max(0, Math.min(estado.duracao, v));
        solicitarRender();
        agendarSalvar();
      });
    }

    if (duracaoInput) {
      var mudarDuracao = function () {
        var v = parseFloat(duracaoInput.value);
        if (isNaN(v)) return;
        v = Math.max(1, Math.min(DURACAO_MAX, v));
        if (v === estado.duracao) return;
        estado.duracao = v;
        aplicarDuracao();
        atualizarLeitura(true);
        agendarSalvar();
      };
      duracaoInput.addEventListener("change", mudarDuracao);
      duracaoInput.addEventListener("input", mudarDuracao);
    }

    if (velocidadeSelect) {
      var achou = false;
      Array.prototype.forEach.call(velocidadeSelect.options, function (op) {
        if (parseFloat(op.value) === estado.velocidade) {
          op.selected = true;
          achou = true;
        }
      });
      if (!achou) estado.velocidade = parseFloat(velocidadeSelect.value) || 60;
      velocidadeSelect.addEventListener("change", function () {
        var v = parseFloat(velocidadeSelect.value);
        if (v > 0) {
          estado.velocidade = v;
          agendarSalvar();
        }
      });
    }

    // Faixas: clique escolhe o tubo, arrasto move o tempo
    var arrasto = null;
    faixasCanvases.forEach(function (canvas) {
      function tempoDoPonteiro(e) {
        var r = canvas.getBoundingClientRect();
        return Math.max(0, Math.min(1, (e.clientX - r.left) / r.width)) * estado.duracao;
      }
      canvas.addEventListener("pointerdown", function (e) {
        if (e.button !== 0) return;
        if (tocando) parar();
        arrasto = { id: e.pointerId, x: e.clientX, moveu: false };
        try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* sem captura */ }
        estado.t = tempoDoPonteiro(e);
        solicitarRender();
      });
      canvas.addEventListener("pointermove", function (e) {
        if (!arrasto || arrasto.id !== e.pointerId) return;
        if (Math.abs(e.clientX - arrasto.x) > 3) arrasto.moveu = true;
        estado.t = tempoDoPonteiro(e);
        solicitarRender();
      });
      function soltar(e) {
        if (!arrasto || arrasto.id !== e.pointerId) return;
        try { canvas.releasePointerCapture(e.pointerId); } catch (err) { /* sem captura */ }
        if (!arrasto.moveu) selecionar(parseInt(canvas.getAttribute("data-tubos-faixa"), 10));
        arrasto = null;
        agendarSalvar();
      }
      canvas.addEventListener("pointerup", soltar);
      canvas.addEventListener("pointercancel", soltar);
      canvas.addEventListener("click", function (e) { e.stopPropagation(); });
    });
    faixasLinhas.forEach(function (linha, i) {
      linha.addEventListener("click", function () { selecionar(i); });
    });

    if (typeof ResizeObserver !== "undefined") {
      var ro = new ResizeObserver(function () {
        sujarFaixas();
        recriarEixo();
        solicitarRender();
      });
      faixasCanvases.forEach(function (c) { ro.observe(c); });
      if (curvaCanvas) ro.observe(curvaCanvas);
    }

    aplicarDuracao();
    atualizarLeitura(true);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", iniciar);
  } else {
    iniciar();
  }
})();
