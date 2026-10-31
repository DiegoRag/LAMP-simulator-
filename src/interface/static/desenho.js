"use strict";

(function () {
  function iniciar() {
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
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", iniciar);
  } else {
    iniciar();
  }
})();
