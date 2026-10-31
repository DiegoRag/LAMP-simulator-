"use strict";

// Fluxo guiado do LampSim: envio automático do arquivo escolhido, bloqueio de duplo envio
// e camada de espera do passo 2 (o envio do passo 2 espera o tempo mínimo do config,
// mesmo com movimento reduzido). Sem dependências; o CSS cuida de toda a animação.
(function () {
  function rotuloDe(form) {
    return form.querySelector("[data-envio-rotulo]");
  }

  function camadaEspera() {
    return document.querySelector("[data-espera-fluxo]");
  }

  function mostrarEspera() {
    var camada = camadaEspera();
    if (!camada) return;
    camada.hidden = false;
    document.documentElement.classList.add("espera-ativa");
    var conteudo = camada.querySelector(".espera-fluxo-conteudo");
    if (conteudo) {
      conteudo.setAttribute("tabindex", "-1");
      conteudo.focus();
    }
  }

  function esconderEspera() {
    var camada = camadaEspera();
    if (camada) camada.hidden = true;
    document.documentElement.classList.remove("espera-ativa");
  }

  function enviar(form) {
    if (form.requestSubmit) form.requestSubmit();
    else form.submit();
  }

  function aoEscolherArquivo(form, entrada) {
    if (!entrada.files || entrada.files.length === 0) return;
    if (form.dataset.enviando) return;
    form.dataset.enviando = "1";
    var rotulo = rotuloDe(form);
    if (rotulo) rotulo.textContent = "Enviando " + entrada.files[0].name + "…";
    if (form.dataset.espera === "fluxo") {
      mostrarEspera();
      var ms = parseInt(form.dataset.esperaMs, 10);
      if (ms > 0) {
        setTimeout(function () {
          enviar(form);
        }, ms);
        return;
      }
    }
    enviar(form);
  }

  function prepararFormulario(form) {
    var rotulo = rotuloDe(form);
    if (rotulo) form.dataset.rotuloOriginal = rotulo.textContent;
    var entrada = form.querySelector('input[type="file"]');
    if (entrada) {
      entrada.addEventListener("change", function () {
        aoEscolherArquivo(form, entrada);
      });
    }
    form.addEventListener("submit", function (e) {
      // um segundo envio só passa se este for o próprio envio já marcado (requestSubmit acima)
      if (form.dataset.enviado) {
        e.preventDefault();
        return;
      }
      if (form.dataset.enviando) form.dataset.enviado = "1";
    });
  }

  function restaurarFormulario(form) {
    delete form.dataset.enviando;
    delete form.dataset.enviado;
    var rotulo = rotuloDe(form);
    if (rotulo && form.dataset.rotuloOriginal !== undefined) {
      rotulo.textContent = form.dataset.rotuloOriginal;
    }
    var entrada = form.querySelector('input[type="file"]');
    if (entrada) entrada.value = "";
  }

  function iniciar() {
    var formularios = document.querySelectorAll("form[data-envio-automatico]");
    Array.prototype.forEach.call(formularios, prepararFormulario);

    window.addEventListener("pageshow", function (e) {
      if (!e.persisted) return;
      Array.prototype.forEach.call(formularios, restaurarFormulario);
      esconderEspera();
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", iniciar);
  } else {
    iniciar();
  }
})();
