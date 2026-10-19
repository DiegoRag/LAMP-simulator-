"use strict";

// Sinaliza para o HTML que o script app.js rodou com sucesso
document.documentElement.classList.add("app-pronto");

// a transição entre páginas pode ser abortada pelo navegador (recarga, aba oculta);
// a falha não importa, então ela não deve aparecer como erro no console
window.addEventListener("pagereveal", function (evento) {
  var vt = evento.viewTransition;
  if (vt) {
    vt.ready.catch(function () {});
    vt.finished.catch(function () {});
  }
});
window.addEventListener("pageswap", function (evento) {
  var vt = evento.viewTransition;
  if (vt) {
    vt.ready.catch(function () {});
    vt.finished.catch(function () {});
  }
});

(function () {
  var REDUZIDO = document.documentElement.classList.contains("movimento-reduzido");

  // Etapas reais do fluxo no servidor
  var ETAPAS = {
    analise: [
      "Lendo o conjunto e o alvo",
      "Localizando as regiões no alvo",
      "Calculando Tm, GC e a ponta 3′",
      "Procurando grampos e dímeros",
      "Conferindo as regras do config",
      "Montando o relatório"
    ],
    desenho: [
      "Gerando candidatos para cada região",
      "Descartando trechos repetidos",
      "Combinando as seis regiões pelas distâncias",
      "Ordenando por Tm e heterodímeros",
      "Encaixando os loops LF e LB",
      "Autoverificando cada proposta"
    ],
    envio: [
      "Recebendo o arquivo",
      "Conferindo o formato",
      "Guardando na memória desta sessão"
    ]
  };

  // Dicas mostradas durante a tela de espera
  var DICAS = [
    "O LAMP copia o DNA a 65 °C constantes, sem termociclador.",
    "FIP e BIP são primers de duas partes: F1c + F2 e B1c + B2.",
    "A Tm é a temperatura em que metade dos primers se solta do alvo.",
    "Quanto mais G e C, mais forte a ligação do primer ao alvo.",
    "Os primers de loop, LF e LB, aceleram a reação.",
    "No ensaio colorimétrico, o tubo passa de rosa a amarelo quando o pH cai.",
    "A referência é a mesma mistura sem a amostra: ela deve continuar rosa.",
    "O gene fosA dá à bactéria resistência à fosfomicina.",
    "Mesma entrada e mesmo config produzem o mesmo resultado, byte a byte.",
    "Nenhuma sequência sai desta máquina: o servidor só escuta em 127.0.0.1."
  ];

  function sgm(x) {
    return 1 / (1 + Math.exp(-x));
  }

  // Curva de amplificação em forma sigmoide (devagar, rápido, platô)
  function curva(u) {
    var a = sgm(-0.42 / 0.12);
    var b = sgm(0.58 / 0.12);
    return (sgm((u - 0.42) / 0.12) - a) / (b - a);
  }

  document.addEventListener("DOMContentLoaded", function () {
    // Alternar animações (movimento reduzido)
    document.querySelectorAll("[data-alternar-movimento]").forEach(function (botao) {
      botao.setAttribute("aria-pressed", REDUZIDO ? "false" : "true");
      botao.addEventListener("click", function () {
        var ligar = botao.getAttribute("aria-pressed") !== "true";
        try {
          localStorage.setItem("lamp-movimento", ligar ? "ligado" : "reduzido");
        } catch (e) {
          // sem armazenamento: a troca vale só para esta carga
        }
        location.reload();
      });
    });

    // Aviso de animações desligadas (rodada 8)
    (function () {
      if (!REDUZIDO) return;

      var semEscolha = false;
      try {
        semEscolha = localStorage.getItem("lamp-movimento") === null;
      } catch (e) {
        semEscolha = false;
      }
      if (!semEscolha) return;

      var sistemaPedeMenos = Boolean(
        window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches
      );
      if (!sistemaPedeMenos) return;

      var barraTopoDireita = document.querySelector(".barra-topo-direita");
      if (!barraTopoDireita) return;

      var dica = document.createElement("div");
      dica.className = "dica-movimento";
      dica.setAttribute("role", "status");

      var p = document.createElement("p");
      p.className = "dica-movimento-texto";
      p.textContent = "As animações estão desligadas porque o sistema pede menos movimento.";
      dica.appendChild(p);

      var acoes = document.createElement("div");
      acoes.className = "dica-movimento-acoes";

      var btnLigar = document.createElement("button");
      btnLigar.type = "button";
      btnLigar.className = "botao botao-primario botao-pequeno";
      btnLigar.setAttribute("data-dica-ligar", "");
      btnLigar.textContent = "Ligar animações";

      var btnFechar = document.createElement("button");
      btnFechar.type = "button";
      btnFechar.className = "botao botao-secundario botao-pequeno";
      btnFechar.setAttribute("data-dica-fechar", "");
      btnFechar.textContent = "Deixar assim";

      function fecharDica() {
        try {
          localStorage.setItem("lamp-movimento", "reduzido");
        } catch (e) {
          // sem armazenamento
        }
        if (typeof dica.remove === "function") {
          dica.remove();
        } else if (dica.parentNode) {
          dica.parentNode.removeChild(dica);
        }
      }

      btnLigar.addEventListener("click", function () {
        try {
          localStorage.setItem("lamp-movimento", "ligado");
        } catch (e) {
          // sem armazenamento
        }
        location.reload();
      });

      btnFechar.addEventListener("click", fecharDica);

      dica.addEventListener("keydown", function (e) {
        if (e.key === "Escape" || e.key === "Esc") {
          fecharDica();
        }
      });

      acoes.appendChild(btnLigar);
      acoes.appendChild(btnFechar);
      dica.appendChild(acoes);

      barraTopoDireita.appendChild(dica);
    })();

    // C1. Entradas em cascata e encaixes
    (function () {
      var elementosRevelar = document.querySelectorAll("[data-revelar]");
      elementosRevelar.forEach(function (el) {
        var passo = parseFloat(el.dataset.revelarPasso) || 0.06;
        var filhos = el.children;
        for (var i = 0; i < filhos.length; i++) {
          var atraso = Math.min(i * passo, 1.2);
          filhos[i].style.setProperty("--d", atraso + "s");
        }
      });

      var observados = document.querySelectorAll("[data-revelar], [data-encaixe], [data-onda]");
      if (REDUZIDO || typeof IntersectionObserver === "undefined") {
        observados.forEach(function (el) {
          el.classList.add("visivel");
        });
      } else {
        var ioRevelar = new IntersectionObserver(
          function (entries) {
            entries.forEach(function (entry) {
              if (entry.isIntersecting) {
                entry.target.classList.add("visivel");
                ioRevelar.unobserve(entry.target);
              }
            });
          },
          { threshold: 0.12, rootMargin: "0px 0px -8% 0px" }
        );
        observados.forEach(function (el) {
          ioRevelar.observe(el);
        });
      }
    })();

    // C2. Números que contam
    (function () {
      if (REDUZIDO) return;

      var regexNumero = /[−-]?\d{1,3}(?:\.\d{3})+(?:,\d+)?|[−-]?\d+(?:,\d+)?/;
      var elementos = document.querySelectorAll("[data-contar]");

      function animarElemento(el) {
        var textoOriginal = el.textContent;
        var match = textoOriginal.match(regexNumero);
        if (!match) return;

        var numStr = match[0];
        var antes = textoOriginal.slice(0, match.index);
        var depois = textoOriginal.slice(match.index + numStr.length);

        var sinal = "";
        var sResto = numStr;
        if (numStr[0] === "−" || numStr[0] === "-") {
          sinal = numStr[0];
          sResto = numStr.slice(1);
        }

        var temMilhar = sResto.indexOf(".") !== -1;
        var partesVirgula = sResto.split(",");
        var casas = partesVirgula.length > 1 ? partesVirgula[1].length : 0;
        var inteiroLimpo = partesVirgula[0].replace(/\./g, "");
        var valorFinal = parseFloat(inteiroLimpo + (casas > 0 ? "." + partesVirgula[1] : ""));
        if (sinal) valorFinal = -valorFinal;

        function formatarNumero(v) {
          var absV = Math.abs(v);
          var fixo = absV.toFixed(casas);
          var partes = fixo.split(".");
          var parteInt = partes[0];
          if (temMilhar) {
            parteInt = parteInt.replace(/\B(?=(\d{3})+(?!\d))/g, ".");
          }
          var parteDec = casas > 0 ? "," + partes[1] : "";
          var prefixoSinal = (v < 0 || (v === 0 && sinal)) ? sinal : "";
          return antes + prefixoSinal + parteInt + parteDec + depois;
        }

        var inicio = null;
        var duracao = 900;

        function quadro(agora) {
          if (!inicio) inicio = agora;
          var decorrido = agora - inicio;
          var prog = Math.min(decorrido / duracao, 1);
          var ease = 1 - Math.pow(1 - prog, 3);
          var vAtual = valorFinal * ease;

          el.textContent = formatarNumero(vAtual);

          if (prog < 1) {
            requestAnimationFrame(quadro);
          } else {
            el.textContent = textoOriginal;
          }
        }

        requestAnimationFrame(quadro);
      }

      if (typeof IntersectionObserver !== "undefined") {
        var ioContar = new IntersectionObserver(
          function (entries) {
            entries.forEach(function (entry) {
              if (entry.isIntersecting) {
                animarElemento(entry.target);
                ioContar.unobserve(entry.target);
              }
            });
          },
          { threshold: 0.12, rootMargin: "0px 0px -8% 0px" }
        );
        elementos.forEach(function (el) {
          ioContar.observe(el);
        });
      } else {
        elementos.forEach(animarElemento);
      }
    })();

    // C3. Hélice da abertura
    (function () {
      var canvasAnatomia = document.querySelector('canvas[data-helice="anatomia"]');
      if (canvasAnatomia && window.LampHelice) {
        var heliceAnatomia = window.LampHelice.criar(canvasAnatomia, { tema: "claro" });
        heliceAnatomia.iniciar();

        if (REDUZIDO) {
          heliceAnatomia.encaixarTudo();
        } else if (typeof IntersectionObserver !== "undefined") {
          var disparou = false;
          var ioAnatomia = new IntersectionObserver(
            function (entries) {
              entries.forEach(function (entry) {
                if (entry.isIntersecting && !disparou) {
                  disparou = true;
                  heliceAnatomia.encaixarEmSequencia(170);
                  ioAnatomia.disconnect();
                }
              });
            },
            { threshold: 0.1 }
          );
          ioAnatomia.observe(canvasAnatomia);
        } else {
          heliceAnatomia.encaixarEmSequencia(170);
        }
      }
    })();

    // Hélice de dados da análise (rodada 5)
    (function () {
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
        var rawRegioes = canvasDados.getAttribute("data-regioes") || canvasDados.dataset.regioes;
        if (!rawRegioes) return;
        var segmentos;
        try {
          segmentos = JSON.parse(rawRegioes);
        } catch (e) {
          return;
        }
        if (!Array.isArray(segmentos)) return;

        var rawComprimento = canvasDados.getAttribute("data-comprimento") || canvasDados.dataset.comprimento;
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
            fantasma: Boolean(s.fantasma)
          });
        }

        var opcoesHelice = {
          tema: "claro",
          regioes: lista
        };
        if (!isNaN(n)) {
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

    // C4. Tela de espera
    var telaEspera = document.getElementById("tela-espera");
    var heliceEspera = null;
    var esperaAberta = false;
    var animacaoProgressoId = null;
    var timerDica = null;
    var timerDestino = null;

    var elTitulo = telaEspera ? telaEspera.querySelector("[data-espera-titulo]") : null;
    var elTag = telaEspera ? telaEspera.querySelector("[data-espera-tag]") : null;
    var elEtapasLista = telaEspera ? telaEspera.querySelector("[data-espera-etapas]") : null;
    var elEtapaTexto = telaEspera ? telaEspera.querySelector("[data-espera-etapa]") : null;
    var elValor = telaEspera ? telaEspera.querySelector("[data-espera-valor]") : null;
    var elDica = telaEspera ? telaEspera.querySelector("[data-espera-dica]") : null;
    var elBarra = telaEspera ? telaEspera.querySelector(".espera-progresso-barra") : null;
    var canvasEspera = telaEspera ? telaEspera.querySelector('canvas[data-helice="espera"]') : null;

    function obterHeliceEspera() {
      if (!heliceEspera && canvasEspera && window.LampHelice) {
        heliceEspera = window.LampHelice.criar(canvasEspera, { tema: "escuro" });
      }
      return heliceEspera;
    }

    function obterTempoMinimo(perfil) {
      var dataset = document.body ? document.body.dataset : {};
      if (perfil === "analise") {
        return parseInt(dataset.esperaAnalise, 10) || 2600;
      }
      if (perfil === "desenho") {
        return parseInt(dataset.esperaDesenho, 10) || 3200;
      }
      return parseInt(dataset.esperaEnvio, 10) || 900;
    }

    function atualizarVisualProgresso(p, listaEtapasTexto) {
      if (!telaEspera) return;

      telaEspera.style.setProperty("--p", p);

      // Mistura linear de rosa (230, 57, 127) para amarelo (245, 184, 0)
      var r = Math.round(230 + (245 - 230) * p);
      var g = Math.round(57 + (184 - 57) * p);
      var b = Math.round(127 + (0 - 127) * p);
      telaEspera.style.setProperty("--cor-progresso", "rgb(" + r + "," + g + "," + b + ")");

      var pct = Math.round(p * 100);
      if (elValor) elValor.textContent = pct + " %";
      if (elBarra) elBarra.setAttribute("aria-valuenow", pct);

      if (elEtapasLista) {
        var itens = elEtapasLista.querySelectorAll(".espera-etapa");
        var n = itens.length;
        var i = Math.min(n - 1, Math.floor(p * n));

        for (var idx = 0; idx < n; idx++) {
          if (p >= 0.999) {
            itens[idx].className = "espera-etapa feita";
          } else if (idx < i) {
            itens[idx].className = "espera-etapa feita";
          } else if (idx === i) {
            itens[idx].className = "espera-etapa atual";
          } else {
            itens[idx].className = "espera-etapa";
          }
        }

        var textoEtapa;
        if (p >= 0.999) {
          textoEtapa = "Pronto.";
        } else {
          textoEtapa = (listaEtapasTexto && listaEtapasTexto[i] ? listaEtapasTexto[i] : "") + "…";
        }
        if (elEtapaTexto && elEtapaTexto.textContent !== textoEtapa) {
          elEtapaTexto.textContent = textoEtapa;
        }
      }
    }

    function fecharEspera() {
      esperaAberta = false;
      if (animacaoProgressoId) {
        cancelAnimationFrame(animacaoProgressoId);
        animacaoProgressoId = null;
      }
      if (timerDica) {
        clearInterval(timerDica);
        timerDica = null;
      }
      if (timerDestino) {
        clearTimeout(timerDestino);
        timerDestino = null;
      }

      if (telaEspera) {
        telaEspera.classList.remove("espera-ativa");
        telaEspera.hidden = true;
        telaEspera.style.setProperty("--p", "0");
        telaEspera.style.removeProperty("--cor-progresso");
      }

      document.body.classList.remove("esperando");
      var barraTopo = document.querySelector(".barra-topo");
      var main = document.querySelector("main");
      if (barraTopo) barraTopo.removeAttribute("inert");
      if (main) main.removeAttribute("inert");

      if (heliceEspera) {
        heliceEspera.parar();
      }
    }

    function abrirEspera(perfil, titulo, tag, callbackDestinoPronto) {
      if (esperaAberta || !telaEspera) return;
      esperaAberta = true;

      var etapas = ETAPAS[perfil] || ETAPAS.envio;
      if (elTitulo) elTitulo.textContent = titulo || "Aguarde";
      if (elTag) elTag.textContent = tag ? tag : "";

      if (elEtapasLista) {
        elEtapasLista.innerHTML = "";
        for (var eIdx = 0; eIdx < etapas.length; eIdx++) {
          var li = document.createElement("li");
          li.className = "espera-etapa" + (eIdx === 0 ? " atual" : "");
          li.textContent = etapas[eIdx];
          elEtapasLista.appendChild(li);
        }
      }

      var indiceDica = Math.floor(Math.random() * DICAS.length);
      if (elDica) elDica.textContent = DICAS[indiceDica];

      telaEspera.hidden = false;
      requestAnimationFrame(function () {
        telaEspera.classList.add("espera-ativa");
      });

      document.body.classList.add("esperando");
      var barraTopo = document.querySelector(".barra-topo");
      var main = document.querySelector("main");
      if (barraTopo) barraTopo.setAttribute("inert", "");
      if (main) main.setAttribute("inert", "");

      var h2 = telaEspera.querySelector("h2");
      if (h2) h2.focus();

      var helice = obterHeliceEspera();
      if (helice) {
        helice.iniciar();
      }

      timerDica = setInterval(function () {
        if (!elDica) return;
        elDica.classList.add("trocando");
        setTimeout(function () {
          indiceDica = (indiceDica + 1) % DICAS.length;
          elDica.textContent = DICAS[indiceDica];
          elDica.classList.remove("trocando");
        }, 350);
      }, 3400);

      var min = obterTempoMinimo(perfil);
      var tempoInicio = performance.now();
      var ultimoTempo = tempoInicio;
      var p = 0;
      var respondeu = false;
      var urlDestino = "";

      function loopProgresso(agora) {
        if (!esperaAberta) return;

        var el = agora - tempoInicio;
        var dt = ultimoTempo ? (agora - ultimoTempo) / 1000 : 0.016;
        ultimoTempo = agora;
        if (dt > 0.05) dt = 0.05;

        var alvo;
        if (respondeu && el >= min) {
          alvo = 1;
        } else if (el <= min) {
          alvo = 0.94 * curva(Math.min(el / Math.max(min, 1), 1));
        } else {
          alvo = 0.94 + 0.05 * (1 - Math.exp(-(el - min) / 4000));
        }

        p += (alvo - p) * (1 - Math.exp(-dt * 8));
        if (alvo === 1 && Math.abs(1 - p) < 0.002) {
          p = 1;
        }

        atualizarVisualProgresso(p, etapas);
        if (helice) {
          helice.definirProgresso(p);
        }

        if (respondeu && el >= min && p >= 0.999) {
          atualizarVisualProgresso(1, etapas);
          if (helice) helice.definirProgresso(1);
          if (elEtapaTexto) elEtapaTexto.textContent = "Pronto.";

          timerDestino = setTimeout(function () {
            if (urlDestino) {
              if (urlDestino === location.pathname + location.search) {
                location.reload();
              } else {
                location.assign(urlDestino);
              }
            }
          }, 320);
          return;
        }

        animacaoProgressoId = requestAnimationFrame(loopProgresso);
      }

      animacaoProgressoId = requestAnimationFrame(loopProgresso);

      if (typeof callbackDestinoPronto === "function") {
        callbackDestinoPronto(function (dest) {
          urlDestino = dest;
          respondeu = true;
        });
      }
    }

    function abrirEsperaCongelada(perfil, titulo, tag, pCongelado) {
      if (!telaEspera) return;
      esperaAberta = true;

      var etapas = ETAPAS[perfil] || ETAPAS.envio;
      if (elTitulo) elTitulo.textContent = titulo;
      if (elTag) elTag.textContent = tag ? tag : "";

      if (elEtapasLista) {
        elEtapasLista.innerHTML = "";
        for (var eIdx = 0; eIdx < etapas.length; eIdx++) {
          var li = document.createElement("li");
          li.className = "espera-etapa";
          li.textContent = etapas[eIdx];
          elEtapasLista.appendChild(li);
        }
      }

      if (elDica) elDica.textContent = DICAS[0];

      telaEspera.hidden = false;
      telaEspera.classList.add("espera-ativa");
      document.body.classList.add("esperando");

      var barraTopo = document.querySelector(".barra-topo");
      var main = document.querySelector("main");
      if (barraTopo) barraTopo.setAttribute("inert", "");
      if (main) main.setAttribute("inert", "");

      var h2 = telaEspera.querySelector("h2");
      if (h2) h2.focus();

      atualizarVisualProgresso(pCongelado, etapas);

      var helice = obterHeliceEspera();
      if (helice) {
        if (pCongelado >= 1) {
          helice.encaixarTudo();
        } else {
          helice.definirProgresso(pCongelado, true);
        }
      }
    }

    function enviarComEspera(form, perfil, titulo, tag) {
      abrirEspera(perfil, titulo, tag, function (marcarPronto) {
        fetch(form.action, {
          method: "POST",
          body: new FormData(form),
          headers: { "X-Requested-With": "lamp-carregamento" },
          credentials: "same-origin"
        })
          .then(function (res) {
            if (!res.ok) {
              throw new Error("Resposta inválida: " + res.status);
            }
            return res.json();
          })
          .then(function (dados) {
            if (dados && dados.destino) {
              marcarPronto(dados.destino);
            } else {
              throw new Error("JSON sem destino");
            }
          })
          .catch(function () {
            fecharEspera();
            form.submit();
          });
      });
    }

    // Intercepta formulários com data-carregar
    document.querySelectorAll("form[data-carregar]").forEach(function (form) {
      form.addEventListener("submit", function (e) {
        if (typeof window.fetch !== "function") return;
        e.preventDefault();

        var perfil = form.dataset.carregar;
        var titulo = form.dataset.titulo || "Aguarde";
        var tag = form.dataset.tag || "";
        if (!tag) {
          var inputArq = form.querySelector('input[type="file"]');
          if (inputArq && inputArq.files && inputArq.files.length > 0) {
            tag = inputArq.files[0].name;
          }
        }
        enviarComEspera(form, perfil, titulo, tag);
      });
    });

    // Fecha a espera ao navegar pelo histórico (bfcache)
    window.addEventListener("pageshow", function (event) {
      if (event.persisted) {
        fecharEspera();
      }
    });

    // Modo de demonstração da tela de espera
    (function () {
      if (typeof window === "undefined" || !window.location) return;
      var params = new URLSearchParams(window.location.search);
      var demoPerfil = params.get("demo-espera");
      if (demoPerfil && (demoPerfil === "analise" || demoPerfil === "desenho" || demoPerfil === "envio")) {
        var titulos = {
          analise: "Análise",
          desenho: "Desenho",
          envio: "Envio"
        };
        var paramCongelar = params.get("congelar");
        if (paramCongelar !== null && !isNaN(parseFloat(paramCongelar))) {
          var pCong = parseFloat(paramCongelar);
          abrirEsperaCongelada(demoPerfil, titulos[demoPerfil], "demonstração", pCong);
        } else {
          abrirEspera(demoPerfil, titulos[demoPerfil], "demonstração");
        }
      }
    })();

    // C5. O que já existia (data-copiar, data-alternar, input[type=file])
    document.querySelectorAll("[data-copiar]").forEach(function (btn) {
      var timeoutId = null;
      btn.addEventListener("click", function () {
        var texto = btn.getAttribute("data-copiar");
        var rotulo = btn.querySelector("span") || btn;

        function definirRotulo(novo) {
          rotulo.textContent = novo;
          if (timeoutId) {
            clearTimeout(timeoutId);
          }
          timeoutId = setTimeout(function () {
            rotulo.textContent = "Copiar";
            timeoutId = null;
          }, 2000);
        }

        if (navigator.clipboard && navigator.clipboard.writeText) {
          navigator.clipboard.writeText(texto).then(
            function () {
              definirRotulo("Copiado");
            },
            function () {
              definirRotulo("Selecione e copie");
            }
          );
        } else {
          definirRotulo("Selecione e copie");
        }
      });
    });

    document.querySelectorAll("[data-alternar]").forEach(function (elemento) {
      elemento.addEventListener("change", function () {
        var destino = elemento.dataset.alternar || elemento.getAttribute("data-alternar");
        if (destino) {
          window.location.href = destino;
        }
      });
    });

    document.querySelectorAll("form input[type=file]").forEach(function (input) {
      input.addEventListener("change", function () {
        if (input.files && input.files.length > 0 && input.form) {
          var indicador = input.form.querySelector("[data-nome-arquivo]");
          if (indicador) {
            indicador.textContent = input.files[0].name;
          }
          if (input.form.dataset.carregar) {
            if (typeof input.form.requestSubmit === "function") {
              input.form.requestSubmit();
            } else {
              var perfil = input.form.dataset.carregar;
              var titulo = input.form.dataset.titulo || "Envio";
              var tag = input.form.dataset.tag || input.files[0].name;
              enviarComEspera(input.form, perfil, titulo, tag);
            }
          } else {
            input.form.submit();
          }
        }
      });
    });
  });
})();
