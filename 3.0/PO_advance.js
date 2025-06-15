const fs = require("fs");

main();

function main(){
    let linesArray = readTxt();

    let constraintCount = countConstraints(linesArray);

    let [xCount, optimizationType] = getXCount(linesArray, constraintCount);

    let inequalityValues = addVariables(linesArray, xCount);

    let [fullMatrix, mainExpressionVector] = fillMatrix(linesArray, xCount, constraintCount);

    console.log('Here is the complete matrix:');
    console.table(fullMatrix);

    let determinant = calculateDeterminant(fullMatrix);

    console.log('Matrix determinant:', determinant);
    console.log('Main expression vector after adding new variables:', mainExpressionVector);

    if (checkPhaseOne(linesArray, mainExpressionVector, inequalityValues, fullMatrix)[0]) {
        phaseOne(fullMatrix, inequalityValues, mainExpressionVector, optimizationType);
    } else {
        const numRows = fullMatrix.length;
        const numCols = fullMatrix[0].length;

        const columnIndices = [...Array(numCols).keys()];

        for (let i = columnIndices.length - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [columnIndices[i], columnIndices[j]] = [columnIndices[j], columnIndices[i]];
        }

        let basicColumns = columnIndices.slice(0, numRows);
        let nonBasicColumns = columnIndices.slice(numRows);

        console.log("Shuffled column indices:", columnIndices);
        console.log("Basic columns:", basicColumns);
        console.log("Non-basic columns:", nonBasicColumns);

        let basicMatrix = fullMatrix.map(row =>
            basicColumns.map(i => row[i])
        );
        console.log(`Basic matrix: `, basicMatrix);

        let nonBasicMatrix = fullMatrix.map(row =>
            nonBasicColumns.map(i => row[i])
        );
        console.log(`Non-basic matrix: `, nonBasicMatrix);

        phaseTwo(fullMatrix, basicMatrix, basicColumns, nonBasicMatrix, nonBasicColumns, inequalityValues, mainExpressionVector, optimizationType);
    }
}

function readTxt(){
    const array = fs.readFileSync("teste.txt", "utf-8")
        .split("\n")
        .map(line => line.trim().replace(/\s+/g, ""))
        .filter(line => line !== "");

    return array;
}

function countConstraints(array){
    let constraintCount = 0;
    for(let i = 1; i < array.length; i++){
        if(array[i].includes(">=") || array[i].includes("<=")){
            constraintCount++;
        }
    }
    console.log("Number of constraint lines: ", constraintCount);

    return constraintCount;
}

function getXCount(array){
    let xCount = 0;
    for(let i = 0; i < array[0].length; i++){
        if(array[0][i] === "x"){
            xCount++;
        }
    }

    let optimizationType;
    if(array[0].toLowerCase().startsWith("max")){
        optimizationType = "max";
        xCount--;
    } else {
        optimizationType = "min";
    }

    console.log("Number of X variables in the main expression: ", xCount);

    return [xCount, optimizationType];
}

function addVariables(array, xCount){
    let artificialIndex = xCount + 1;
    let inequalityValues = [];

    for (let i = 1; i < array.length; i++) {
        if (array[i].includes(">=")) {
            let [left, right] = array[i].split(">=");
            inequalityValues.push(eval(right));
            array[i] = `${left}-x${artificialIndex}>=${right}`;
            artificialIndex++;
        } else if (array[i].includes("<=")) {
            let [left, right] = array[i].split("<=");
            inequalityValues.push(eval(right));
            array[i] = `${left}+x${artificialIndex}<=${right}`;
            artificialIndex++;
        } else if (array[i].includes("=") && !array[i].includes(">=") && !array[i].includes("<=")) {
            let [left, right] = array[i].split("=");
            inequalityValues.push(eval(right));
        }
    }

    console.log(`Values after inequality (vector b):`, inequalityValues);
    console.log("Updated array:", array);

    return inequalityValues;
}

function fillMatrix(array, xCount, constraintCount){
    let fullMatrix = Array(array.length - 1).fill(null).map(() => Array(xCount + constraintCount).fill(0));
    let regex = /([+-]?\d*\.?\d*)\*?x(\d+)/g;

    let mainExpressionVector = [];

    for (let i = 0; i < array.length; i++) {
        let expr = array[i];

        for (let match of expr.matchAll(regex)) {
            let coef = match[1];
            let xIndex = parseInt(match[2]) - 1;

            if (coef === "" || coef === "+") {
                coef = 1;
            } else if (coef === "-") {
                coef = -1;
            }

            if(i === 0){
                mainExpressionVector.push(parseFloat(coef));
            } else {
                fullMatrix[i - 1][xIndex] = parseFloat(coef);
            }
        }
    }

    for(let i = 0; i < constraintCount; i++){
        mainExpressionVector.push(0);
    }

    return [fullMatrix, mainExpressionVector];
}

function calculateDeterminant(matrix) {
    const n = matrix.length;

    if (n === 1) return matrix[0][0];

    let det = 0;
    const row = 0;

    for (let col = 0; col < n; col++) {
        det += ((-1) ** (row + col)) * matrix[row][col] * calculateDeterminant(subMatrix(matrix, row, col));
    }
    return det;
}

function subMatrix(matrix, row, col) {
    return matrix
        .filter((_, i) => i !== row)
        .map(currentRow => currentRow.filter((_, j) => j !== col));
}

function createSquareMatrix(fullMatrix){
    let n = fullMatrix.length;
    let squareMatrix = Array.from({ length: n }, () => Array(n).fill(0));

    for(let i = 0; i < n; i++){
        for(let j = 0; j < n; j++){
            squareMatrix[i][j] = fullMatrix[i][j];
        }
    }
    return squareMatrix;
}

function createIdentityMatrix(squareMatrix){
    let identity = Array.from({ length: squareMatrix.length }, () => Array(squareMatrix.length).fill(0));

    for (let i = 0; i < identity.length; i++) {
        for (let j = 0; j < identity.length; j++){
            if(i === j){
                identity[i][i] = 1;
            }
        }
    }
    return identity;
}

function createInverseMatrix(squareMatrix, identity){
    let augmentedMatrix = squareMatrix.map((row, i) => row.concat(identity[i]));
    console.log('Augmented matrix = squareMatrix + identity');
    console.table(augmentedMatrix);

    for (let i = 0; i < squareMatrix.length; i++) {
        if (augmentedMatrix[i][i] === 0) {
            let swapped = false;
            for (let k = i + 1; k < squareMatrix.length; k++) {
                if (augmentedMatrix[k][i] !== 0) {
                    [augmentedMatrix[i], augmentedMatrix[k]] = [augmentedMatrix[k], augmentedMatrix[i]];
                    swapped = true;
                    break;
                }
            }
            if (!swapped) {
                console.log("Matrix is not invertible");
                return null;
            }
        }

        let pivot = augmentedMatrix[i][i];
        for (let j = 0; j < 2 * squareMatrix.length; j++) {
            augmentedMatrix[i][j] /= pivot;
        }

        for (let k = 0; k < squareMatrix.length; k++) {
            if (k !== i) {
                let factor = augmentedMatrix[k][i];
                for (let j = 0; j < 2 * squareMatrix.length; j++) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }
    }

    let inverse = augmentedMatrix.map(row => row.slice(squareMatrix.length));
    return inverse;
}

function multiplyMatrix(matrix1, matrix2){
    let finalMatrix = Array(matrix1.length).fill(0).map(() => Array(matrix2[0].length).fill(0));
    if(matrix1[0].length === matrix2.length){
        for(let i = 0; i < matrix1.length; i++){
            for(let j = 0; j < matrix2[0].length; j++){
                let sum = 0;
                for(let k = 0; k < matrix1[0].length; k++){
                    sum += matrix1[i][k] * matrix2[k][j];
                }
                finalMatrix[i][j] = sum;
            }
        }
    } else {
        console.log("Incompatible matrix sizes for multiplication");
    }
    return finalMatrix;
}

function createBasicMatrix(matrix) {
    const numRows = matrix.length;
    const numCols = matrix[0].length;

    const columnIndices = [...Array(numCols).keys()];
    let basicColumns = columnIndices.slice(0, numRows);
    let nonBasicColumns = columnIndices.slice(numRows);

    let basicMatrix = matrix.map(row =>
        basicColumns.map(i => row[i])
    );
    console.log(basicMatrix);

    let nonBasicMatrix = matrix.map(row =>
        nonBasicColumns.map(i => row[i])
    );
    console.log(nonBasicMatrix);

    let attempt = 0;
    while (calculateDeterminant(basicMatrix) === 0 && attempt < 100) {
        console.log("Initial basic matrix has determinant 0. Trying another basis...");

        for (let i = columnIndices.length - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [columnIndices[i], columnIndices[j]] = [columnIndices[j], columnIndices[i]];
        }

        console.log("Shuffled column indices: ", columnIndices);

        basicColumns = columnIndices.slice(0, numRows);
        basicMatrix = matrix.map(row =>
            basicColumns.map(i => row[i])
        );

        nonBasicColumns = columnIndices.slice(numRows);
        nonBasicMatrix = matrix.map(row =>
            nonBasicColumns.map(i => row[i])
        );

        attempt++;
    }

    if (calculateDeterminant(basicMatrix) === 0) {
        console.log("No basic matrix found with non-zero determinant");
        return null;
    }

    return [basicMatrix, nonBasicMatrix, basicColumns, nonBasicColumns];
}

function checkPhaseOne(array, mainExpressionVector, inequalityValues, fullMatrix){
    if (array[0].toLowerCase().startsWith("max")) {
        for(let i = 0; i < mainExpressionVector.length; i++){
            mainExpressionVector[i] *= -1;
        }
    }

    for(let i = 0; i < inequalityValues.length; i++){
        if(inequalityValues[i] < 0){
            inequalityValues[i] *= -1;
            for(let j = 0; j < fullMatrix[i].length; j++){
                fullMatrix[i][j] = fullMatrix[i][j] * -1;
            }
        }
    }

    for (let i = 1; i < array.length; i++) {
        if (array[i].includes(">=") || array[i].includes(">")) {
            console.log("Going to Phase I");
            return [true, mainExpressionVector, fullMatrix];
        }

        if (array[i].includes("=") && !array[i].includes(">=") && !array[i].includes("<=")) {
            console.log("Going to Phase I");
            return [true, mainExpressionVector, fullMatrix];
        }
    }

    console.log("Going to Phase II");
    return [false, mainExpressionVector, fullMatrix];
}

function phaseTwo(fullMatrix, basicMatrix, basicColumns, nonBasicMatrix, nonBasicColumns, inequalityValues, mainExpressionVector, optimizationType){
    let iteration = 1;
    while(iteration < 100){
        console.log(`Phase II, iteration: ${iteration}`);
        
        console.table(fullMatrix);
        console.log(basicColumns)
        console.log(nonBasicColumns)
        console.table(basicMatrix);
        console.table(nonBasicMatrix);

        let inverseBasic = createInverseMatrix(basicMatrix, createIdentityMatrix(basicMatrix));
        console.log("Inverse of basic matrix:");
        console.table(inverseBasic);

        let vectorB = inequalityValues.map(i => [i]);
        let xBasic = multiplyMatrix(inverseBasic, vectorB);
        console.log("Basic solution vector:");
        console.table(xBasic);

        let costBasic = [basicColumns.map(i => mainExpressionVector[i])];
        console.log("Basic cost vector:", costBasic);
        let yt = multiplyMatrix(costBasic, inverseBasic);
        console.log("yt vector:");
        console.table(yt);

        let costNonBasic = [nonBasicColumns.map(i => mainExpressionVector[i])];
        console.log("Non-basic cost vector:", costNonBasic);

        let aNj = Array(fullMatrix.length).fill().map(() => 
            Array(nonBasicColumns.length).fill(0)
        );
        for(let i = 0; i < fullMatrix.length; i++) {
            for(let j = 0; j < nonBasicColumns.length; j++) {
                aNj[i][j] = fullMatrix[i][nonBasicColumns[j]];
            }
        }
        console.log("aNj matrix:");
        console.table(aNj);

        let mult = multiplyMatrix(yt, aNj);
        console.log("Multiplication yt * aNj:");
        console.table(mult);

        let relativeCost = [];
        for(let i = 0; i < costNonBasic[0].length; i++){
            relativeCost.push(costNonBasic[0][i] - mult[0][i]);
        }
        console.log("Relative cost:");
        console.table(relativeCost);

        let enteringIndex = relativeCost.indexOf(Math.min(...relativeCost));
        let enteringValue = Math.min(...relativeCost);
        console.log("Entering variable index:");
        console.log(enteringIndex);
        console.log("Minimum relative cost value:");
        console.log(enteringValue);

        if(enteringValue >= 0){
            console.log(`Optimal solution found at iteration: ${iteration}`);

            let optimalValue = 0;
            for(let i = 0; i < costBasic[0].length; i++){
                optimalValue += costBasic[0][i] * xBasic[i][0];
            }
            if(optimizationType === "max"){
                optimalValue *= -1;
            }

            let solutionVector = Array(mainExpressionVector.length).fill(0);
            for(let i = 0; i < xBasic.length; i++){
                const varIndex = basicColumns[i];
                solutionVector[varIndex] = xBasic[i][0];
            }

            console.log(solutionVector);
            console.log(`Optimal value of the objective function (z): ${optimalValue}`);

            return [
                optimalValue,
                solutionVector,
                iteration
            ];
        }

        let aNk = aNj.map(row => [row[enteringIndex]]);
        console.log("aNk column vector:");
        console.table(aNk);

        let direction = multiplyMatrix(inverseBasic, aNk);
        console.log("Simplex direction (y):");
        console.table(direction);

        if(direction.every(element => element <= 0)){
            console.log("Unbounded problem. No finite optimal solution.");
            return null;
        }

        let epsilon = Infinity;
        let leavingIndex = -1;
        for(let i = 0; i < direction.length; i++){
            if(direction[i][0] > 0){
                let ratio = xBasic[i][0] / direction[i][0];
                if(ratio < epsilon){
                    epsilon = ratio;
                    leavingIndex = i;
                    console.log(`Epsilon (ε̂): ${epsilon}`);
                    console.log(`Index of variable leaving the basis: ${leavingIndex}`);
                }
            }
        }

        let enteringVar = nonBasicColumns[enteringIndex];
        let leavingVar = basicColumns[leavingIndex];
        basicColumns[leavingIndex] = enteringVar;
        nonBasicColumns[enteringIndex] = leavingVar;
        console.log(basicColumns);
        console.log(nonBasicColumns);

        basicMatrix = fullMatrix.map(row =>
            basicColumns.map(i => row[i])
        );
        console.table(basicMatrix);

        nonBasicMatrix = fullMatrix.map(row =>
            nonBasicColumns.map(i => row[i])
        );
        console.table(nonBasicMatrix);

        console.log(`Variable entering the basis: x${enteringVar + 1}`);
        console.log(`Variable leaving the basis: x${leavingVar + 1}`);

        iteration++;
    }
}

function phaseOne(fullMatrix, inequalityValues, mainExpressionVector, originalBasicColumns, optimizationType){
    let m = fullMatrix.length;
    let n = fullMatrix[0].length;

    let matrixWithArtificial = fullMatrix.map((row, i) => {
        let artificial = Array(m).fill(0);
        artificial[i] = 1;
        return [...row, ...artificial];
    });

    let artificialExpression = Array(n).fill(0).concat(Array(m).fill(1));

    let basicColumns = [];
    let nonBasicColumns = [];
    for(let i = 0; i < n + m; i++){
        if(i >= n){
            basicColumns.push(i);
        } else {
            nonBasicColumns.push(i);
        }
    }

    let basicMatrix = matrixWithArtificial.map(row => basicColumns.map(i => row[i]));
    let nonBasicMatrix = matrixWithArtificial.map(row => nonBasicColumns.map(i => row[i]));

    console.table(basicMatrix);
    console.table(nonBasicMatrix);

    let resultPhaseOne = phaseTwo(
        matrixWithArtificial,
        basicMatrix,
        basicColumns,
        nonBasicMatrix,
        nonBasicColumns,
        inequalityValues,
        artificialExpression,
        optimizationType
    );

    if(resultPhaseOne === null || resultPhaseOne[0] > 1e-6){
        console.log("Infeasible problem. No initial feasible basic solution.");
        return null;
    }

    let matrixWithoutArtificial = matrixWithArtificial.map(row => row.slice(0, n));

    let newBasicColumns = basicColumns.filter(c => c < n);
    let newNonBasicColumns = nonBasicColumns.filter(c => c < n);
    while(newBasicColumns.length < m){
        let candidate = newNonBasicColumns.shift();
        newBasicColumns.push(candidate);
    }

    let newBasicMatrix = matrixWithoutArtificial.map(row => newBasicColumns.map(i => row[i]));
    let newNonBasicMatrix = matrixWithoutArtificial.map(row => newNonBasicColumns.map(i => row[i]));

    return phaseTwo(
        matrixWithoutArtificial,
        newBasicMatrix,
        newBasicColumns,
        newNonBasicMatrix,
        newNonBasicColumns,
        inequalityValues,
        mainExpressionVector,
        optimizationType
    );
}
