const fileSystem = require("fs");

main();

function main(){
    let inputLines = readTextFile();

    let totalConstraints = countAllConstraints(inputLines);

    let [variableCount, problemType] = getVariableCount(inputLines, totalConstraints);

    let constraintValues = insertAdditionalVariables(inputLines, variableCount);

    let [completeMatrix, objectiveVector] = buildCompleteMatrix(inputLines, variableCount, totalConstraints);

    console.log('Here is the complete matrix:');
    console.table(completeMatrix);

    let matrixDet = computeMatrixDeterminant(completeMatrix);

    console.log('Matrix determinant:', matrixDet);
    console.log('Main expression vector after adding new variables:', objectiveVector);

    if (verifyFirstPhaseNeeded(inputLines, objectiveVector, constraintValues, completeMatrix)[0]) {
        runFirstPhase(completeMatrix, constraintValues, objectiveVector, problemType);
    } else {
        const rowCount = completeMatrix.length;
        const columnCount = completeMatrix[0].length;

        const allColumns = Array.from({length: columnCount}, (_, i) => i);
        const randomSeed = Date.now() % 1000;
        allColumns.sort((a, b) => {
            const hash1 = (a * 2654435761 ^ randomSeed) % 1000;
            const hash2 = (b * 2654435761 ^ randomSeed) % 1000;
            return hash1 - hash2;
        });

        let basisColumns = allColumns.slice(0, rowCount);
        let nonBasisColumns = allColumns.slice(rowCount);

        console.log("Shuffled column indices:", allColumns);
        console.log("Basic columns:", basisColumns);
        console.log("Non-basic columns:", nonBasisColumns);

        let basisSubmatrix = completeMatrix.map(row =>
            basisColumns.map(i => row[i])
        );
        console.log(`Basic matrix: `, basisSubmatrix);

        let nonBasisSubmatrix = completeMatrix.map(row =>
            nonBasisColumns.map(i => row[i])
        );
        console.log(`Non-basic matrix: `, nonBasisSubmatrix);

        runSecondPhase(completeMatrix, basisSubmatrix, basisColumns, nonBasisSubmatrix, nonBasisColumns, constraintValues, objectiveVector, problemType);
    }
}

function readTextFile(){
    const data = fileSystem.readFileSync("TXT.txt", "utf-8")
        .split("\n")
        .map(line => line.trim().replace(/\s+/g, ""))
        .filter(line => line !== "");

    return data;
}

function countAllConstraints(dataArray){
    let constraintTotal = 0;
    for(let i = 1; i < dataArray.length; i++){
        if(dataArray[i].includes(">=") || dataArray[i].includes("<=")){
            constraintTotal++;
        }
    }
    console.log("Number of constraint lines: ", constraintTotal);

    return constraintTotal;
}

function getVariableCount(dataArray){
    let varCount = 0;
    for(let i = 0; i < dataArray[0].length; i++){
        if(dataArray[0][i] === "x"){
            varCount++;
        }
    }

    let optimizationGoal;
    if(dataArray[0].toLowerCase().startsWith("max")){
        optimizationGoal = "max";
        varCount--;
    } else {
        optimizationGoal = "min";
    }

    console.log("Number of X variables in the main expression: ", varCount);

    return [varCount, optimizationGoal];
}

function insertAdditionalVariables(dataArray, varCount){
    let newVarIndex = varCount + 1;
    let constraintValues = [];

    for (let i = 1; i < dataArray.length; i++) {
        if (dataArray[i].includes(">=")) {
            let [leftPart, rightPart] = dataArray[i].split(">=");
            constraintValues.push(eval(rightPart));
            dataArray[i] = `${leftPart}-x${newVarIndex}>=${rightPart}`;
            newVarIndex++;
        } else if (dataArray[i].includes("<=")) {
            let [leftPart, rightPart] = dataArray[i].split("<=");
            constraintValues.push(eval(rightPart));
            dataArray[i] = `${leftPart}+x${newVarIndex}<=${rightPart}`;
            newVarIndex++;
        } else if (dataArray[i].includes("=") && !dataArray[i].includes(">=") && !dataArray[i].includes("<=")) {
            let [leftPart, rightPart] = dataArray[i].split("=");
            constraintValues.push(eval(rightPart));
        }
    }

    console.log(`Values after inequality (vector b):`, constraintValues);
    console.log("Updated array:", dataArray);

    return constraintValues;
}

function buildCompleteMatrix(dataArray, varCount, constraintTotal){
    let finalMatrix = Array(dataArray.length - 1).fill(null).map(() => Array(varCount + constraintTotal).fill(0));
    let pattern = /([+-]?\d*\.?\d*)\*?x(\d+)/g;

    let objectiveCoefficients = [];

    for (let i = 0; i < dataArray.length; i++) {
        let expression = dataArray[i];

        for (let match of expression.matchAll(pattern)) {
            let coefficient = match[1];
            let varPosition = parseInt(match[2]) - 1;

            if (coefficient === "" || coefficient === "+") {
                coefficient = 1;
            } else if (coefficient === "-") {
                coefficient = -1;
            }

            if(i === 0){
                objectiveCoefficients.push(parseFloat(coefficient));
            } else {
                finalMatrix[i - 1][varPosition] = parseFloat(coefficient);
            }
        }
    }

    for(let i = 0; i < constraintTotal; i++){
        objectiveCoefficients.push(0);
    }

    return [finalMatrix, objectiveCoefficients];
}

function computeMatrixDeterminant(inputMatrix) {
    const size = inputMatrix.length;

    if (size === 1) return inputMatrix[0][0];

    let determinant = 0;
    const currentRow = 0;

    for (let col = 0; col < size; col++) {
        determinant += ((-1) ** (currentRow + col)) * inputMatrix[currentRow][col] * computeMatrixDeterminant(getSubmatrix(inputMatrix, currentRow, col));
    }
    return determinant;
}

function getSubmatrix(matrix, excludeRow, excludeCol) {
    return matrix
        .filter((_, rowIndex) => rowIndex !== excludeRow)
        .map(currentRow => currentRow.filter((_, colIndex) => colIndex !== excludeCol));
}

function makeSquareMatrix(fullMatrix){
    let size = fullMatrix.length;
    let squareMatrix = Array.from({ length: size }, () => Array(size).fill(0));

    for(let i = 0; i < size; i++){
        for(let j = 0; j < size; j++){
            squareMatrix[i][j] = fullMatrix[i][j];
        }
    }
    return squareMatrix;
}

function generateIdentityMatrix(squareMatrix){
    let identityMatrix = Array.from({ length: squareMatrix.length }, () => Array(squareMatrix.length).fill(0));

    for (let i = 0; i < identityMatrix.length; i++) {
        for (let j = 0; j < identityMatrix.length; j++){
            if(i === j){
                identityMatrix[i][i] = 1;
            }
        }
    }
    return identityMatrix;
}

function computeInverseMatrix(squareMatrix, identityMatrix){
    let combinedMatrix = squareMatrix.map((row, i) => row.concat(identityMatrix[i]));
    console.log('Augmented matrix = squareMatrix + identity');
    console.table(combinedMatrix);

    for (let i = 0; i < squareMatrix.length; i++) {
        if (combinedMatrix[i][i] === 0) {
            let foundSwap = false;
            for (let k = i + 1; k < squareMatrix.length; k++) {
                if (combinedMatrix[k][i] !== 0) {
                    [combinedMatrix[i], combinedMatrix[k]] = [combinedMatrix[k], combinedMatrix[i]];
                    foundSwap = true;
                    break;
                }
            }
            if (!foundSwap) {
                console.log("Matrix is not invertible");
                return null;
            }
        }

        let pivotValue = combinedMatrix[i][i];
        for (let j = 0; j < 2 * squareMatrix.length; j++) {
            combinedMatrix[i][j] /= pivotValue;
        }

        for (let k = 0; k < squareMatrix.length; k++) {
            if (k !== i) {
                let factor = combinedMatrix[k][i];
                for (let j = 0; j < 2 * squareMatrix.length; j++) {
                    combinedMatrix[k][j] -= factor * combinedMatrix[i][j];
                }
            }
        }
    }

    let inverseMatrix = combinedMatrix.map(row => row.slice(squareMatrix.length));
    return inverseMatrix;
}

function matrixProduct(matrixA, matrixB){
    let resultMatrix = Array(matrixA.length).fill(0).map(() => Array(matrixB[0].length).fill(0));
    if(matrixA[0].length === matrixB.length){
        for(let i = 0; i < matrixA.length; i++){
            for(let j = 0; j < matrixB[0].length; j++){
                let sum = 0;
                for(let k = 0; k < matrixA[0].length; k++){
                    sum += matrixA[i][k] * matrixB[k][j];
                }
                resultMatrix[i][j] = sum;
            }
        }
    } else {
        console.log("Incompatible matrix sizes for multiplication");
    }
    return resultMatrix;
}

function initializeBasis(matrix) {
    const rowNum = matrix.length;
    const colNum = matrix[0].length;

    const allIndices = [...Array(colNum).keys()];
    let basisIndices = allIndices.slice(0, rowNum);
    let nonBasisIndices = allIndices.slice(rowNum);

    let basisMat = matrix.map(row =>
        basisIndices.map(i => row[i])
    );
    console.log(basisMat);

    let nonBasisMat = matrix.map(row =>
        nonBasisIndices.map(i => row[i])
    );
    console.log(nonBasisMat);

    let attempts = 0;
    while (computeMatrixDeterminant(basisMat) === 0 && attempts < 100) {
        console.log("Initial basic matrix has determinant 0. Trying another basis...");

        for (let i = allIndices.length - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [allIndices[i], allIndices[j]] = [allIndices[j], allIndices[i]];
        }

        console.log("Shuffled column indices: ", allIndices);

        basisIndices = allIndices.slice(0, rowNum);
        basisMat = matrix.map(row =>
            basisIndices.map(i => row[i])
        );

        nonBasisIndices = allIndices.slice(rowNum);
        nonBasisMat = matrix.map(row =>
            nonBasisIndices.map(i => row[i])
        );

        attempts++;
    }

    if (computeMatrixDeterminant(basisMat) === 0) {
        console.log("No basic matrix found with non-zero determinant");
        return null;
    }

    return [basisMat, nonBasisMat, basisIndices, nonBasisIndices];
}

function verifyFirstPhaseNeeded(dataArray, objectiveVector, constraintValues, fullMatrix){
    if (dataArray[0].toLowerCase().startsWith("max")) {
        for(let i = 0; i < objectiveVector.length; i++){
            objectiveVector[i] *= -1;
        }
    }

    for(let i = 0; i < constraintValues.length; i++){
        if(constraintValues[i] < 0){
            constraintValues[i] *= -1;
            for(let j = 0; j < fullMatrix[i].length; j++){
                fullMatrix[i][j] = fullMatrix[i][j] * -1;
            }
        }
    }

    for (let i = 1; i < dataArray.length; i++) {
        if (dataArray[i].includes(">=") || dataArray[i].includes(">")) {
            console.log("Going to Phase I");
            return [true, objectiveVector, fullMatrix];
        }

        if (dataArray[i].includes("=") && !dataArray[i].includes(">=") && !dataArray[i].includes("<=")) {
            console.log("Going to Phase I");
            return [true, objectiveVector, fullMatrix];
        }
    }

    console.log("Going to Phase II");
    return [false, objectiveVector, fullMatrix];
}

function runSecondPhase(fullMatrix, basisMat, basisIndices, nonBasisMat, nonBasisIndices, constraintValues, objectiveVector, problemType){
    let step = 1;
    while(step < 100){
        console.log(`Phase II, iteration: ${step}`);
        
        console.table(fullMatrix);
        console.log(basisIndices)
        console.log(nonBasisIndices)
        console.table(basisMat);
        console.table(nonBasisMat);

        let inverseBasis = computeInverseMatrix(basisMat, generateIdentityMatrix(basisMat));
        console.log("Inverse of basic matrix:");
        console.table(inverseBasis);

        let bVector = constraintValues.map(i => [i]);
        let xBasis = matrixProduct(inverseBasis, bVector);
        console.log("Basic solution vector:");
        console.table(xBasis);

        let cBasis = [basisIndices.map(i => objectiveVector[i])];
        console.log("Basic cost vector:", cBasis);
        let yTranspose = matrixProduct(cBasis, inverseBasis);
        console.log("yt vector:");
        console.table(yTranspose);

        let cNonBasis = [nonBasisIndices.map(i => objectiveVector[i])];
        console.log("Non-basic cost vector:", cNonBasis);

        let aNonBasis = Array(fullMatrix.length).fill().map(() => 
            Array(nonBasisIndices.length).fill(0)
        );
        for(let i = 0; i < fullMatrix.length; i++) {
            for(let j = 0; j < nonBasisIndices.length; j++) {
                aNonBasis[i][j] = fullMatrix[i][nonBasisIndices[j]];
            }
        }
        console.log("aNj matrix:");
        console.table(aNonBasis);

        let product = matrixProduct(yTranspose, aNonBasis);
        console.log("Multiplication yt * aNj:");
        console.table(product);

        let reducedCost = [];
        for(let i = 0; i < cNonBasis[0].length; i++){
            reducedCost.push(cNonBasis[0][i] - product[0][i]);
        }
        console.log("Relative cost:");
        console.table(reducedCost);

        let enteringVarIndex = reducedCost.indexOf(Math.min(...reducedCost));
        let minReducedCost = Math.min(...reducedCost);
        console.log("Entering variable index:");
        console.log(enteringVarIndex);
        console.log("Minimum relative cost value:");
        console.log(minReducedCost);

        if(minReducedCost >= 0){
            console.log(`Optimal solution found at iteration: ${step}`);

            let optimalSolution = 0;
            for(let i = 0; i < cBasis[0].length; i++){
                optimalSolution += cBasis[0][i] * xBasis[i][0];
            }
            if(problemType === "max"){
                optimalSolution *= -1;
            }

            let solution = Array(objectiveVector.length).fill(0);
            for(let i = 0; i < xBasis.length; i++){
                const varIndex = basisIndices[i];
                solution[varIndex] = xBasis[i][0];
            }

            console.log(solution);
            console.log(`Optimal value of the objective function (z): ${optimalSolution}`);

            return [
                optimalSolution,
                solution,
                step
            ];
        }

        let enteringColumn = aNonBasis.map(row => [row[enteringVarIndex]]);
        console.log("aNk column vector:");
        console.table(enteringColumn);

        let directionVector = matrixProduct(inverseBasis, enteringColumn);
        console.log("Simplex direction (y):");
        console.table(directionVector);

        if(directionVector.every(element => element <= 0)){
            console.log("Unbounded problem. No finite optimal solution.");
            return null;
        }

        let minRatio = Infinity;
        let leavingVarIndex = -1;
        for(let i = 0; i < directionVector.length; i++){
            if(directionVector[i][0] > 0){
                let currentRatio = xBasis[i][0] / directionVector[i][0];
                if(currentRatio < minRatio){
                    minRatio = currentRatio;
                    leavingVarIndex = i;
                    console.log(`Epsilon (ε̂): ${minRatio}`);
                    console.log(`Index of variable leaving the basis: ${leavingVarIndex}`);
                }
            }
        }

        let enteringVariable = nonBasisIndices[enteringVarIndex];
        let leavingVariable = basisIndices[leavingVarIndex];
        basisIndices[leavingVarIndex] = enteringVariable;
        nonBasisIndices[enteringVarIndex] = leavingVariable;
        console.log(basisIndices);
        console.log(nonBasisIndices);

        basisMat = fullMatrix.map(row =>
            basisIndices.map(i => row[i])
        );
        console.table(basisMat);

        nonBasisMat = fullMatrix.map(row =>
            nonBasisIndices.map(i => row[i])
        );
        console.table(nonBasisMat);

        console.log(`Variable entering the basis: x${enteringVariable + 1}`);
        console.log(`Variable leaving the basis: x${leavingVariable + 1}`);

        step++;
    }
}

function runFirstPhase(fullMatrix, constraintValues, objectiveVector, originalBasisIndices, problemType){
    let rows = fullMatrix.length;
    let cols = fullMatrix[0].length;

    let extendedMatrix = fullMatrix.map((row, i) => {
        let artificialVars = Array(rows).fill(0);
        artificialVars[i] = 1;
        return [...row, ...artificialVars];
    });

    let phaseOneObjective = Array(cols).fill(0).concat(Array(rows).fill(1));

    let basisCols = [];
    let nonBasisCols = [];
    for(let i = 0; i < cols + rows; i++){
        if(i >= cols){
            basisCols.push(i);
        } else {
            nonBasisCols.push(i);
        }
    }

    let basisMatrix = extendedMatrix.map(row => basisCols.map(i => row[i]));
    let nonBasisMatrix = extendedMatrix.map(row => nonBasisCols.map(i => row[i]));

    console.table(basisMatrix);
    console.table(nonBasisMatrix);

    let phaseOneResult = runSecondPhase(
        extendedMatrix,
        basisMatrix,
        basisCols,
        nonBasisMatrix,
        nonBasisCols,
        constraintValues,
        phaseOneObjective,
        problemType
    );

    if(phaseOneResult === null || phaseOneResult[0] > 1e-6){
        console.log("Infeasible problem. No initial feasible basic solution.");
        return null;
    }

    let reducedMatrix = extendedMatrix.map(row => row.slice(0, cols));

    let newBasisCols = basisCols.filter(c => c < cols);
    let newNonBasisCols = nonBasisCols.filter(c => c < cols);
    while(newBasisCols.length < rows){
        let candidate = newNonBasisCols.shift();
        newBasisCols.push(candidate);
    }

    let newBasisMat = reducedMatrix.map(row => newBasisCols.map(i => row[i]));
    let newNonBasisMat = reducedMatrix.map(row => newNonBasisCols.map(i => row[i]));

    return runSecondPhase(
        reducedMatrix,
        newBasisMat,
        newBasisCols,
        newNonBasisMat,
        newNonBasisCols,
        constraintValues,
        objectiveVector,
        problemType
    );
}
