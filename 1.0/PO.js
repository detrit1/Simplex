const fs = require("fs");

main();

function readTxtFile() {
    const lines = fs.readFileSync("TXT.txt", "utf-8")
        .split("\n") 
        .map(line => line.trim())
        .filter(line => line !== ""); 

    for (let i = 0; i < lines.length; i++) {
        if (!lines[i].includes("=")) {
            lines.splice(i, 1);
        }
    }

    if (lines[0].toLowerCase().startsWith("max")) {
        console.log("É max");
    } else if (lines[0].toLowerCase().startsWith("min")) {
        console.log("É min");
    } else {
        console.log("Você deve especificar se é 'min' ou 'max' no início do arquivo");
    }

    parseConstraints(lines);
    return lines;
}

function parseConstraints(lines) {
    let constraintCount = 0;
    for (let i = 1; i < lines.length; i++) {
        if (lines[i].includes(">=") || lines[i].includes("<=")) {
            constraintCount++;
        }
    }
    console.log("Quantidade de restrições:", constraintCount);

    countVariables(lines, constraintCount);
}

function countVariables(lines, constraintCount) {
    let variableCount = 0;
    for (let i = 0; i < lines[0].length; i++) {
        if (lines[0][i] === "x") {
            variableCount++;
        }
    }
    if (lines[0].toLowerCase().startsWith("max")) {
        variableCount--;
    }
    console.log("Quantidade de variáveis na função objetivo:", variableCount);
    addSlackVariables(lines, variableCount, constraintCount);
}

function addSlackVariables(lines, variableCount, constraintCount) {
    let auxVarIndex = variableCount + 1;
    let rhsValues = [];

    for (let i = 1; i < lines.length; i++) {
        if (lines[i].includes(">=")) {
            let [lhs, rhs] = lines[i].split(">=");
            rhsValues.push(rhs);
            lines[i] = `${lhs}-x${auxVarIndex}>=${rhs}`;
            auxVarIndex++;
        } else if (lines[i].includes("<=")) {
            let [lhs, rhs] = lines[i].split("<=");
            rhsValues.push(rhs);
            lines[i] = `${lhs}+x${auxVarIndex}<=${rhs}`;
            auxVarIndex++;
        }
    }

    console.log("Valores à direita das desigualdades:", rhsValues);
    console.log("Array atualizado:", lines);
    buildMatrix(lines, variableCount, constraintCount);
}

function buildMatrix(lines, variableCount, constraintCount) {
    let fullMatrix = Array(lines.length - 1).fill(null).map(() => Array(variableCount + constraintCount).fill(0));

    let regex = /([+-]?\d*\.?\d*)\*?x(\d+)/g;
    let objectiveVector = [];

    for (let i = 0; i < lines.length; i++) {
        let expression = lines[i];
        for (let match of expression.matchAll(regex)) {
            let coef = match[1];
            let xIndex = parseInt(match[2]) - 1;

            if (coef === "" || coef === "+") coef = 1;
            else if (coef === "-") coef = -1;

            if (i === 0) {
                objectiveVector.push(parseFloat(coef));
            } else {
                fullMatrix[i - 1][xIndex] = parseFloat(coef);
            }
        }
    }

    for (let i = 0; i < constraintCount; i++) {
        objectiveVector.push(0);
    }

    console.log("Matriz completa:");
    console.table(fullMatrix);
    console.log("Determinante da matriz:", determinant(fullMatrix));
    console.log("Vetor objetivo:", objectiveVector);

    let squareMatrix = createSquareMatrix(fullMatrix);
    console.table(squareMatrix);

    let identity = createIdentityMatrix(squareMatrix);
    console.table(identity);

    let inverseMatrix = invertMatrix(squareMatrix, identity);
    console.table(inverseMatrix);

    generateBasicMatrix(fullMatrix);

    let multipliedMatrix = multiplyMatrices(squareMatrix, identity);
    console.table(multipliedMatrix);
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

function createIdentityMatrix(squareMatrix) {
    let identity = Array.from({ length: squareMatrix.length }, () => Array(squareMatrix.length).fill(0));
    for (let i = 0; i < identity.length; i++) {
        identity[i][i] = 1;
    }
    console.log("Matriz identidade:");
    return identity;
}

function invertMatrix(squareMatrix, identity) {
    let augmented = squareMatrix.map((row, i) => row.concat(identity[i]));
    console.log("Matriz aumentada:");
    console.table(augmented);

    for (let i = 0; i < squareMatrix.length; i++) {
        if (augmented[i][i] === 0) {
            let swapped = false;
            for (let k = i + 1; k < squareMatrix.length; k++) {
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

        let factor = augmented[i][i];
        for (let j = 0; j < 2 * squareMatrix.length; j++) {
            augmented[i][j] /= factor;
        }

        for (let k = 0; k < squareMatrix.length; k++) {
            if (k !== i) {
                let factor2 = augmented[k][i];
                for (let j = 0; j < 2 * squareMatrix.length; j++) {
                    augmented[k][j] -= factor2 * augmented[i][j];
                }
            }
        }
    }

    let inverse = augmented.map(row => row.slice(squareMatrix.length));
    console.log("Matriz inversa:");
    return inverse;
}

function generateBasicMatrix(matrix) {
    let columnIndexes = [...Array(matrix[0].length).keys()];
    for (let i = columnIndexes.length - 1; i > 0; i--) {
        const seedrandom = require('seedrandom');
        const rng = seedrandom('seed');
        const j = Math.floor(rng() * (i + 1));
        [columnIndexes[i], columnIndexes[j]] = [columnIndexes[j], columnIndexes[i]];
    }
    console.log("Colunas embaralhadas:", columnIndexes);

    const basicColumns = columnIndexes.slice(0, matrix.length);
    const nonBasicColumns = columnIndexes.slice(matrix.length);

    const basicMatrix = matrix.map(row => basicColumns.map(idx => row[idx]));
    const nonBasicMatrix = matrix.map(row => nonBasicColumns.map(idx => row[idx]));

    console.log("Matriz Básica (B):");
    console.table(basicMatrix);
    console.log("Matriz Não-Básica (N):");
    console.table(nonBasicMatrix);
}

function multiplyMatrices(mat1, mat2) {
    let result = Array(mat1.length).fill(null).map(() => Array(mat2[0].length).fill(0));
    if (mat1[0].length === mat2.length) {
        for (let i = 0; i < mat1.length; i++) {
            for (let j = 0; j < mat2[0].length; j++) {
                for (let k = 0; k < mat1[0].length; k++) {
                    result[i][j] += mat1[i][k] * mat2[k][j];
                }
            }
        }
    } else {
        console.log("Não é possível multiplicar as matrizes.");
    }
    console.log("Matriz resultante da multiplicação:");
    return result;
}

function createBasicMatrix(matrix){
    let columnIndex = [...Array(matrix[0].length).keys()];

    console.log('\níndice organizado de colunas, para da matriz básica e não básica: ', columnIndex);
    console.log("matriz completa:");
    console.table(matrix);

    const columnsToBasic = columnIndexes.slice(0, matrix.length);
    const columnsToNotBasic = columnIndex.slice(matrix.length);

    const notBasicMatrix = matrix.map(line => columnsToBasic.map(index => line[index]));
    const basicMatrix = matrix.map(line => columnsToNotBasic.map(index => line[index]));

    return [basicMatrix, notBasicMatrix, columnsToBasic, columnsToNotBasic];
}

function phase1verify(array, mainExpressionVector, inequalityValues, completeMatrix){
    if(array[0].toLowerCase().startsWith("max")){
        for(let i = 0; i < inequalityValues; i++){
            inequalityValues[i] *= -1;
            for(let j = 0; j < completeMatrix[i].length; j++){
                completeMatrix[i][j] = completeMatrix[i][j] * -1;
            }

        }
    }    
    for(let i = 1; i < array.length; i++){
        if(array[i].includes(">=") || array[i].includes(">")){
            console.log("procede para a fase 1.");
            return[true, mainExpressionVector, completeMatrix]
        }

        if(array[i].includes("=") && array[i].includes(">=") && array[i].includes("<=")){
            console.log("procede para a fase 1.");
            return[true, mainExpressionVector, completeMatrix];
        }
    }
    console.log("procede para a fase 2.");
    return [false, mainExpressionVector, completeMatrix];
}


function pt1(){
    console.log("fase one");
}

function pt2(){
    console.log("fase tu")
}

function main(){
    let array = readTxtFile();
    let lineCounter = parseConstraints(array);
    let [xCounter, optimizationType] = countVariables(array, lineCounter);
    let inequalityValues = addSlackVariables(array, xCounter);
    let [completeMatrix, mainExpressionVector] = buildMatrix(array, xCounter, lineCounter);
    console.log('matriz completa: ')
    console.table(completeMatrix);
    let det = determinant(completeMatrix);
    let [basicMatrix, notBasicMatrix, basicColumns, columnsToNotBasic] = createBasicMatrix(completeMatrix);
    /*if(phase1verify(array, mainExpressionVector, inequalityValues, completeMatrix)[0]){

    }
    else{
        pt2(completeMatrix, inequalityValues, mainExpressionVector, optimizationType);
    }*/
   phase1verify(array, mainExpressionVector, inequalityValues, completeMatrix)[0];
}