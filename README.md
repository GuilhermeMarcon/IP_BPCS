# PDI_BPCS
Guilherme dos Santos Marcon NUSP 9293564

Esteganografia BPCS

Ocultamento de imagens em imagens pelo método BPCS (Bit-Plane Complexity Segmentation)

  O objetivo desse projeto é implementar intuitivamente o método BPCS de esteganografia, que nesse caso servirá para ocultar uma imagem em outra e a recuperar. O método consiste em ocultar pixels da imagem alvo em blocos de um plano de bit cujos bits em si possuem comportamento ruidoso, ele se aproveita na característica da visão humana de se concentrar no reconhecimento de padrões e formas.

Todas as imagens utilizadas estão no diretório "imagens" e elas foram retiradas e modificadas do site https://www.pexels.com/public-domain-images/, sendo transformadas para o formato png e reduzidas o tamanho para testes mais simples. Os links de cada imagem individual estão salvas no arquivo ImagensLinks.txt.

Os principais métodos utilizados são: ler, salvar e manipular as imagens utilizando as bibliotecas imageio e numpy do python; transformar a imagem de Pure Binary Code para Canonical Gray Code e vice-versa; checar se um bloco de um plano de bit é considerado complexo.
