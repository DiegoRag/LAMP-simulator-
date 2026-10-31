"use strict";

(function () {
  var sistemaPedeMenos = Boolean(
    window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches
  );
  if (sistemaPedeMenos && !document.documentElement.classList.contains("movimento-reduzido")) {
    document.documentElement.classList.add("movimento-reduzido");
  }

  var REDUZIDO = document.documentElement.classList.contains("movimento-reduzido");

  var canvasesDados = document.querySelectorAll('canvas[data-helice="dados"]');
  if (!canvasesDados.length || !window.LampHelice) return;

  function familiaPorRegiao(regiao) {
    var r = String(regiao || "").trim().toUpperCase();
    if (r === "F3" || r === "B3C" || r === "B3") return "externo";
    if (r === "F2") return "fip";
    if (r === "F1" || r === "F1C") return "fipClaro";
    if (r === "B1C" || r === "B1") return "bipClaro";
    if (r === "B2C" || r === "B2") return "bip";
    if (r === "LF" || r === "LB") return "loop";
    return "outro";
  }

  canvasesDados.forEach(function (canvasDados) {
    var rawRegioes =
      canvasDados.getAttribute("data-regioes") ||
      (canvasDados.dataset ? canvasDados.dataset.regioes : null);
    if (!rawRegioes) return;
    var segmentos;
    try {
      segmentos = JSON.parse(rawRegioes);
    } catch (e) {
      return;
    }
    if (!Array.isArray(segmentos)) return;

    var rawComprimento =
      canvasDados.getAttribute("data-comprimento") ||
      (canvasDados.dataset ? canvasDados.dataset.comprimento : null);
    var n = parseInt(rawComprimento, 10);

    var lista = [];
    for (var i = 0; i < segmentos.length; i++) {
      var s = segmentos[i];
      if (!s) continue;
      var esq = parseFloat(s.esquerda);
      var larg = parseFloat(s.largura);
      var ini = esq / 100;
      var fim = (esq + larg) / 100;
      lista.push({
        nome: s.nome != null ? s.nome : (s.regiao != null ? s.regiao : ""),
        ini: ini,
        fim: fim,
        familia: familiaPorRegiao(s.regiao || s.nome),
        fantasma: Boolean(s.fantasma),
      });
    }

    var opcoesHelice = {
      tema: "claro",
      regioes: lista,
    };
    if (!isNaN(n) && n > 0) {
      opcoesHelice.comprimento = n;
    }

    var heliceDados = window.LampHelice.criar(canvasDados, opcoesHelice);
    heliceDados.iniciar();

    if (REDUZIDO) {
      heliceDados.encaixarTudo();
    } else if (typeof IntersectionObserver !== "undefined") {
      var disparou = false;
      var ioDados = new IntersectionObserver(
        function (entries) {
          entries.forEach(function (entry) {
            if (entry.isIntersecting && !disparou) {
              disparou = true;
              heliceDados.encaixarEmSequencia(150);
              ioDados.disconnect();
            }
          });
        },
        { threshold: 0.1 }
      );
      ioDados.observe(canvasDados);
    } else {
      heliceDados.encaixarEmSequencia(150);
    }
  });
})();
