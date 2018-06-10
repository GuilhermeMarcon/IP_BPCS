# PDI_BPCS
Guilherme dos Santos Marcon NUSP 9293564

Esteganografia BPCS

Ocultamento de imagens em imagens pelo método BPCS (Bit-Plane Complexity Segmentation)

  O objetivo desse projeto é implementar intuitivamente o método BPCS de esteganografia, que nesse caso servirá para ocultar uma imagem em outra e a recuperar. O método consiste em ocultar pixels da imagem alvo em blocos de um plano de bit cujos bits em si possuem comportamento ruidoso, ele se aproveita na característica da visão humana de se concentrar no reconhecimento de padrões e formas.

Todas as imagens utilizadas estão no diretório "imagens" e elas foram ou retiradas do site https://www.pexels.com/public-domain-images/ ou modificadas a partir das imagens retiradas do mesmo site.

Os principais métodos utilizados serão:

Ler, salvar e manipular as imagens utilizando as bibliotecas imageio e numpy do python.

Transformar a imagem de Pure Binary Code para Canonical Gray Code e vice-versa.

Checar se um bloco de um plano de bit é considerado complexo.
