const fs = require("fs");
const seedrandom = require('seedrandom');

main();

function readTxtFile() {
    const lines = fs.readFileSync("TXT.txt", "utf-8")
        .split("\n") 
        .map(line => line.trim()) 
        .filter(line => line !== ""); 

    for (let i = lines.length - 1; i >= 0; i--) {
        if (!lines[i].includes("=")) {
            lines.splice(i, 1);
        }
    }

    if (lines.length === 0) {
        throw new Error("O arquivo está vazio ou não contém linhas válidas");
    }

    if (lines[0].toLowerCase().startsWith("max")) {
        console.log("É max");
    } else if (lines[0].toLowerCase().startsWith("min")) {
        console.log("É min");
    } else {
        throw new Error("Você deve especificar se é 'min' ou 'max' no início do arquivo");
    }

    return lines;
}

function parseConstraints(lines) {
    let constraintCount = 0;
    for (let i = 1; i < lines.length; i++) {
        if (lines[i].includes(">=") || lines[i].includes("<=") || lines[i].includes("=")) {
            constraintCount++;
        }
    }
    console.log("Quantidade de restrições:", constraintCount);
    return constraintCount;
}

function countVariables(lines, constraintCount) {
    let variableCount = 0;
    const objectiveFunction = lines[0];
    
    const regex = /x\d+/g;
    const matches = objectiveFunction.match(regex);
    
    if (matches) {
        const maxIndex = Math.max(...matches.map(m => parseInt(m.substring(1))));
        variableCount = maxIndex;
    }

    console.log("Quantidade de variáveis na função objetivo:", variableCount);
    return { variableCount, optimizationType: lines[0].toLowerCase().startsWith("max") ? "max" : "min" };
}

function addSlackVariables(lines, variableCount, constraintCount) {
    let auxVarIndex = variableCount + 1;
    let rhsValues = [];

    for (let i = 1; i < lines.length; i++) {
        if (lines[i].includes(">=")) {
            let [lhs, rhs] = lines[i].split(">=");
            rhsValues.push(parseFloat(rhs.trim()));
            lines[i] = `${lhs.trim()}-x${auxVarIndex}>=${rhs.trim()}`;
            auxVarIndex++;
        } else if (lines[i].includes("<=")) {
            let [lhs, rhs] = lines[i].split("<=");
            rhsValues.push(parseFloat(rhs.trim()));
            lines[i] = `${lhs.trim()}+x${auxVarIndex}<=${rhs.trim()}`;
            auxVarIndex++;
        } else if (lines[i].includes("=")) {
            let [lhs, rhs] = lines[i].split("=");
            rhsValues.push(parseFloat(rhs.trim()));
        }
    }

    console.log("Valores à direita das desigualdades:", rhsValues);
    console.log("Array atualizado:", lines);
    return rhsValues;
}

function buildMatrix(lines, variableCount, constraintCount) {
    let fullMatrix = Array(constraintCount).fill(null).map(() => Array(variableCount + constraintCount).fill(0));
    let objectiveVector = Array(variableCount + constraintCount).fill(0);

    const regex = /([+-]?\d*\.?\d+)\*?x(\d+)/g;

    const objectiveFunction = lines[0];
    let match;
    while ((match = regex.exec(objectiveFunction)) !== null) {
        let coef = match[1];
        let xIndex = parseInt(match[2]) - 1;

        if (coef === "" || coef === "+") coef = 1;
        else if (coef === "-") coef = -1;

        objectiveVector[xIndex] = parseFloat(coef);
    }

    for (let i = 1; i < lines.length; i++) {
        const constraint = lines[i];
        regex.lastIndex = 0; 
        while ((match = regex.exec(constraint)) !== null) {
            let coef = match[1];
            let xIndex = parseInt(match[2]) - 1;

            if (coef === "" || coef === "+") coef = 1;
            else if (coef === "-") coef = -1;

            fullMatrix[i-1][xIndex] = parseFloat(coef);
        }
    }

    console.log("Matriz completa:");
    console.table(fullMatrix);
    console.log("Vetor objetivo:", objectiveVector);

    return { completeMatrix: fullMatrix, mainExpressionVector: objectiveVector };
}

function subMatrix(matrix, row, col) {
    return matrix
        .filter((_, i) => i !== row)
        .map(r => r.filter((_, j) => j !== col));
}

function determinant(matrix) {
    const n = matrix.length;
    if (n === 1) return matrix[0][0];
    if (n === 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    let det = 0;
    for (let col = 0; col < n; col++) {
        det += ((-1) ** col) * matrix[0][col] * determinant(subMatrix(matrix, 0, col));
    }
    return det;
}

function createSquareMatrix(matrix) {
    let n = matrix.length;
    let square = Array.from({ length: n }, () => Array(n).fill(0));

    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            square[i][j] = matrix[i][j];
        }
    }
    console.log("Matriz quadrada:");
    return square;
}

function createIdentityMatrix(size) {
    let identity = Array.from({ length: size }, () => Array(size).fill(0));
    for (let i = 0; i < identity.length; i++) {
        identity[i][i] = 1;
    }
    console.log("Matriz identidade:");
    return identity;
}

function invertMatrix(squareMatrix, identity) {
    const n = squareMatrix.length;
    let augmented = squareMatrix.map((row, i) => row.concat(identity[i]));
    
    console.log("Matriz aumentada:");
    console.table(augmented);

    for (let i = 0; i < n; i++) {
        if (augmented[i][i] === 0) {
            let swapped = false;
            for (let k = i + 1; k < n; k++) {
                if (augmented[k][i] !== 0) {
                    [augmented[i], augmented[k]] = [augmented[k], augmented[i]];
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                console.log("A matriz não é invertível.");
                return null;
            }
        }

        let pivot = augmented[i][i];
        for (let j = 0; j < 2 * n; j++) {
            augmented[i][j] /= pivot;
        }

        for (let k = 0; k < n; k++) {
            if (k !== i && augmented[k][i] !== 0) {
                let factor = augmented[k][i];
                for (let j = 0; j < 2 * n; j++) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
            }
        }
    }

    let inverse = augmented.map(row => row.slice(n));
    console.log("Matriz inversa:");
    return inverse;
}

function generateBasicMatrix(matrix) {
    let columnIndexes = [...Array(matrix[0].length).keys()];
    
    const rng = seedrandom('seed');
    for (let i = columnIndexes.length - 1; i > 0; i--) {
        const j = Math.floor(rng() * (i + 1));
        [columnIndexes[i], columnIndexes[j]] = [columnIndexes[j], columnIndexes[i]];
    }
    console.log("Colunas embaralhadas:", columnIndexes);

    const columnsToBasic = columnIndexes.slice(0, matrix.length);
    const columnsToNotBasic = columnIndexes.slice(matrix.length);

    const basicMatrix = matrix.map(row => columnsToBasic.map(idx => row[idx]));
    const nonBasicMatrix = matrix.map(row => columnsToNotBasic.map(idx => row[idx]));

    console.log("Matriz Básica (B):");
    console.table(basicMatrix);
    console.log("Matriz Não-Básica (N):");
    console.table(nonBasicMatrix);
    
    return { basicMatrix, nonBasicMatrix, columnsToBasic, columnsToNotBasic };
}

function multiplyMatrices(mat1, mat2) {
    if (mat1[0].length !== mat2.length) {
        throw new Error("Número de colunas da primeira matriz deve ser igual ao número de linhas da segunda matriz");
    }

    let result = Array(mat1.length).fill(null).map(() => Array(mat2[0].length).fill(0));
    
    for (let i = 0; i < mat1.length; i++) {
        for (let j = 0; j < mat2[0].length; j++) {
            for (let k = 0; k < mat1[0].length; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }
    
    console.log("Matriz resultante da multiplicação:");
    return result;
}

function phase1Verify(lines, mainExpressionVector, inequalityValues, completeMatrix) {
    const optimizationType = lines[0].toLowerCase().startsWith("max") ? "max" : "min";
    
    if (optimizationType === "max") {
        for (let i = 0; i < inequalityValues.length; i++) {
            inequalityValues[i] *= -1;
            for (let j = 0; j < completeMatrix[i].length; j++) {
                completeMatrix[i][j] *= -1;
            }
        }
    }
    
    for (let i = 1; i < lines.length; i++) {
        if (lines[i].includes(">=") || lines[i].includes(">")) {
            console.log("Procede para a fase 1.");
            return { needsPhase1: true, mainExpressionVector, completeMatrix };
        }

        if (lines[i].includes("=") && !lines[i].includes(">=") && !lines[i].includes("<=")) {
            console.log("Procede para a fase 1.");
            return { needsPhase1: true, mainExpressionVector, completeMatrix };
        }
    }
    
    console.log("Procede para a fase 2.");
    return { needsPhase1: false, mainExpressionVector, completeMatrix };
}

function phase1() {
    console.log("Executando fase 1");
}

function phase2(completeMatrix, basicMatrix, columnsToBasic, nonBasicMatrix, columnsToNotBasic, inequalityValues, mainExpressionVector, optimizationType){
    let iteracao = 1;
    while(iteracao < 100){
        console.log(`faseII, iteração: ${iteracao}`);
        
        console.table(completeMatrix);
        console.log(columnsToBasic)
        console.log(columnsToNotBasic)
        console.table(basicMatrix);
        console.table(nonBasicMatrix);
    
        let inversaBasica = invertMatrix(basicMatrix, createIdentityMatrix(basicMatrix));
        console.log("inversa da basica ai o")
        console.table(inversaBasica);
        let vetorB = inequalityValues.map(i => [i])
        let xBasico = multiplyMatrices(inversaBasica, vetorB);
        console.table("x basico:");
        console.table(xBasico);
    
        let custoBasico = [columnsToBasic.map(i => mainExpressionVector[i])];
        console.log("custo basico: ", custoBasico);
        let yt = multiplyMatrices(custoBasico, inversaBasica);
        console.log("o tal do yt, nao faço a minima ideia do q isso significa")
        console.table(yt)
    
        let custoNaoBasico = [columnsToNotBasic.map(i => mainExpressionVector[i])];
        console.log("custo nao basico: ", custoNaoBasico);

        let aNj = Array(completeMatrix.length).fill().map(() => 
            Array(columnsToNotBasic.length).fill(0)
        );
        for(let i = 0; i < completeMatrix.length; i++) {
            for(let j = 0; j < columnsToNotBasic.length; j++) {
                aNj[i][j] = completeMatrix[i][columnsToNotBasic[j]];
            }
        }
        console.log("valores muito locos de aNj aqui");
        console.table(aNj);
        
        let multiplicacao = multiplyMatrices(yt, aNj);
        console.log("multiplicação muito loca de yt * aNj")
        console.table(multiplicacao)
        let custoRelativo = [];
        for(let i = 0; i < custoNaoBasico[0].length; i++){
            custoRelativo.push(custoNaoBasico[0][i] - multiplicacao[0][i]);
        }
        console.table("custoRelativo:");
        console.table(custoRelativo);

        let indiceVariavelEntrada = custoRelativo.indexOf(Math.min(...custoRelativo));
        let variavelEntrada = (Math.min(...custoRelativo));
        console.log("indice da variavel de entrada:")
        console.log(indiceVariavelEntrada)
        console.log("valor minimo das variaveis:")
        console.log(variavelEntrada)

        if(variavelEntrada >= 0){
            console.log(`Solução ótima encontrada na iteração: ${iteracao}`);

            let valorOtimo = 0;
            for(let i = 0; i < custoBasico[0].length; i++){
                valorOtimo += custoBasico[0][i] * xBasico[i][0];
            }
            if(optimizationType === "max"){
                valorOtimo *= -1;
            }

            let vetorSolucao = Array(mainExpressionVector.length).fill(0);
            for(let i = 0; i < xBasico.length; i++){
                const indiceVariavel = columnsToBasic[i];
                vetorSolucao[indiceVariavel] = xBasico[i][0];
            }
            
            console.log(vetorSolucao)
            console.log(`Valor ótimo da função objetivo (z): ${valorOtimo}`);

            return [
                valorOtimo,
                vetorSolucao,
                iteracao
            ];

        }

        let aNk = aNj.map(linha => [linha[indiceVariavelEntrada]]);
        console.table(aNk);
        
        let y = multiplyMatrices(inversaBasica, aNk);
        console.table("calculo da direção simplex (y):");
        console.table(y);

        if(y.every(elemento => elemento <= 0)){
            console.log(`problema não tem solucão ótima finita`);
            console.log(y.map(elemento => elemento<=0));
            return null;
        }
        let epsilon = Infinity;
        let indiceSaida = -1
        for(let i = 0; i < y.length; i++){
            if(y[i][0] > 0){
                let razao = xBasico[i][0]/y[i][0];
                console.log(razao)
                if(razao < epsilon){
                    epsilon = razao;
                    indiceSaida = i;
                    console.log(`Epsilon (ε̂): ${epsilon}`);
                    console.log(`Índice da variável que sai da base: ${indiceSaida}`);
                }
            } else{
                console.log(`nao deu pra dividir nessa iteração ${i+1} por causa de: ${y[i][0]}`);
            }
        }

        
        let entrando = columnsToNotBasic[indiceVariavelEntrada];
        let saindo = columnsToBasic[indiceSaida];
        columnsToBasic[indiceSaida] = entrando;
        columnsToNotBasic[indiceVariavelEntrada] = saindo;
        console.log(columnsToBasic)
        console.log(columnsToNotBasic)

        basicMatrix = completeMatrix.map(linha => 
            columnsToBasic.map(i => linha[i])
        );
        console.table(basicMatrix)

        nonBasicMatrix = completeMatrix.map(linha => 
            columnsToNotBasic.map(i => linha[i])
        );
        console.table(nonBasicMatrix)

        console.log(`Variável que entra na base: x${entrando + 1}`);
        console.log(`Variável que sai da base: x${saindo + 1}`);

        iteracao++;
    }
}


function main() {
    try {
        let lines = readTxtFile();
        let constraintCount = parseConstraints(lines);
        let { variableCount, optimizationType } = countVariables(lines, constraintCount);
        let inequalityValues = addSlackVariables(lines, variableCount, constraintCount);
        let { completeMatrix, mainExpressionVector } = buildMatrix(lines, variableCount, constraintCount);
        let { basicMatrix, nonBasicMatrix, columnsToBasic, columnsToNotBasic } = generateBasicMatrix(completeMatrix);

        
        console.log('Matriz completa:');
        console.table(completeMatrix);
        
        let det = determinant(completeMatrix);
        console.log("Determinante:", det);

        console.log('Matriz completa:');
        console.table(completeMatrix);
        
        console.log('Colunas básicas:', columnsToBasic);
        console.log('Colunas não básicas:', columnsToNotBasic);
        
        console.log('Matriz básica:');
        console.table(basicMatrix);
        
        console.log('Matriz não básica:');
        console.table(nonBasicMatrix);
        
        
        let { needsPhase1 } = phase1Verify(
            lines, 
            mainExpressionVector, 
            inequalityValues, 
            completeMatrix
        );
        
        if (needsPhase1) {
            phase1();
        } else {
            phase2(completeMatrix, basicMatrix, columnsToBasic, nonBasicMatrix, columnsToNotBasic, inequalityValues, mainExpressionVector, optimizationType);
        }
    } catch (error) {
        console.error("Erro:", error.message);
    }
}

