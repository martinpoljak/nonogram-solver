// Nonogram Solver
// alpha version, 20080610
//
// (c) 2008 Martin Kozak (martinkozak@martinkozak.net)
//
// All rights granted also to:
//  Faculty of Nuclear Sciences and Physical Engineering
//  Czech Technical University in Prague
//
// Distributed under the terms of MIT License.
//

program NonogramSolver;

{$APPTYPE CONSOLE}

uses
  SysUtils, Classes, Math;

  ////////// TYPES

  // GENERAL

  type T1DByteArray = Array of Byte;
  type TByteSet = Set of Byte;

  type T1DIntegerArray = Array of Integer;
  type T2DIntegerArray = Array of Array of Integer;
  type T3DIntegerArray = Array of Array of Array of Integer;

  // CONFIGURATION

  type TApplicationParameters = record
    input: String;
    output: String;
    interactive: Boolean;
    verbose: Byte;
  end;

  type TCellRange = Array[(fromCell, toCell)] of Integer;
  type TColumnDataBlockState = (blockInitial, blockStarted, blockFinished);
  type TColumnDataBlock = record
    length: Integer;
    range: TCellRange;
    result: TCellRange;
    state: TColumnDataBlockState;
    final: Integer;
    splitted: Boolean;
  end;
  type TColumnData = Array of TColumnDataBlock;

  type TColumnDataWrapper = record
    length: Integer;
    data: TColumnData;
  end;

  type TDimensionData = Array of TColumnDataWrapper;
  type TData = Array of Array of TColumnDataWrapper;

  type TInitialConfiguration = record
    size: T1DIntegerArray;
    matrix: TData;
    dimensions: Integer;
    final: Integer;
  end;

  // CELLS

  type TBlockList = TByteSet;
  type TCellMatrix = (cellUnknown, cellEmpty, cellFull, cellFinal);

  type TCellMatrixWrapper = record
    state: TCellMatrix;
    blocks: T1DIntegerArray;
    possibleBlocks: Array of TBlockList;
  end;

  type TCellMatrixWrapperPointer = ^TCellMatrixWrapper;

  type TColumnMatrix = Array of TCellMatrixWrapperPointer;
  type TDimensionMatrix = Array of TColumnMatrix;

  type T1DCellMatrixWrapperArray = Array of TCellMatrixWrapper;

  type TDimensionMatrixWrapper = record
    cells: T1DCellMatrixWrapperArray;
  end;



  // ###############################################################

  ////////// FUNCTIONS

  function LoadConfiguration(var input: String) : TInitialConfiguration;
    var fileData : TStringList;
    var i, j : Integer;
    var rowCounter: Integer;
    var currentString : String;
    var currentColumn : ^TColumnDataWrapper;
    var declarationDone : Boolean;
    var finalResult : TInitialConfiguration;
    var parcialStrings : TStringList;
    var currentSize, maxSize : Integer;
    var dimensionCount : Integer;
    var currentDimension : Integer;
    var loadedIndex : Array of Integer;
    var stringLength: Integer;
  begin

    // Reads file
    fileData := TStringList.Create;
    fileData.LoadFromFile(input);
    declarationDone := false;
    rowCounter := 0;

    // Process file
    for i := 0 to (fileData.Count - 1) do begin

      currentString := Trim(fileData[i]);
      stringLength := Length(currentString);

      if (stringLength > 0) and not (currentString[1] = '#') then begin
        if declarationDone = false then
          begin
            declarationDone := true;

            parcialStrings := TStringList.Create;
            parcialStrings.Clear;
            parcialStrings.Delimiter := ' ';
            parcialStrings.DelimitedText := currentString;

            dimensionCount := parcialStrings.Count;

            SetLength(finalResult.size, dimensionCount);
            SetLength(loadedIndex, dimensionCount);
            loadedIndex[0] := 0;
            currentDimension := 0;
            maxSize := 0;

            for j := 0 to (dimensionCount - 1) do begin
              currentSize := StrToInt(parcialStrings[j]);
              finalResult.size[j] := currentSize;
              if currentSize > maxSize then
                maxSize := currentSize;
            end;

//O            SetLength(finalResult.matrix, dimensionCount, maxSize, Ceil(maxSize / 2) + 1);
            SetLength(finalResult.matrix, dimensionCount, maxSize);
            finalResult.dimensions := dimensionCount;
          end

        else if stringLength > 0 then
          begin
            if loadedIndex[currentDimension] = finalResult.size[currentDimension] then
              begin
                Inc(currentDimension);
                loadedIndex[currentDimension] := 0;
                rowCounter := 0;
              end
            else
              Inc(loadedIndex[currentDimension]);

            parcialStrings.Clear;
            parcialStrings.Delimiter := ' ';
            parcialStrings.DelimitedText := currentString;

            SetLength(finalResult.matrix[currentDimension][rowCounter].data, parcialStrings.Count + 1);
            currentColumn := @finalResult.matrix[currentDimension][rowCounter];

            if not (parcialStrings[0] = '0') then
              begin
                for j := 0 to (parcialStrings.Count - 1) do begin
                  currentColumn.data[j].length := StrToInt(parcialStrings[parcialStrings.Count - j - 1]);
                  currentColumn.data[j].state := blockInitial;
                  currentColumn.data[j].final := 0;
                  currentColumn.data[j].result[fromCell] := -1;
                  currentColumn.data[j].result[toCell] := -1;
                end;
                currentColumn.length := parcialStrings.Count;
                currentColumn.data[j].length := 0;
              end
            else
              begin
                currentColumn.length := 0;
                currentColumn.data[0].length := 0;
              end;

            Inc(rowCounter);

          end;
      end;
    end;

    finalResult.final := 0;
    LoadConfiguration := finalResult;

  end;

  // CELL MATRIX TO CHARACTER ######################################

  function CellMatrixToChr(cellData: TCellMatrix): String;
    var finalResult: String;
  begin
    case cellData of
      cellEmpty: finalResult := '.';
      cellFull: finalResult := 'X';
      cellUnknown: finalResult := '?';
      cellFinal: finalResult := '#';
    end;
    CellMatrixToChr := finalResult;
  end;

  // GETS MATRIX CELL ##########################################
  function GetMatrixCell(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; dimension: Integer; column: Integer; cell: Integer): TCellMatrixWrapperPointer;
    var i: Integer;
    var columnSize: Integer;
  begin
    if not (dimension = 0) then begin
      i := column;
      column := cell;
      cell := i;
    end;

    columnSize := configuration.size[1];
    i := (columnSize * column) + cell;
    GetMatrixCell := @matrix.cells[i];
  end;

  // GETS MATRIX COLUMMN #######################################
  function GetMatrixColumn(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; dimension: Integer; column: Integer): TColumnMatrix;
    var i: Integer;
    var step, last: Integer;
    var columnSize: Integer;
    var finalResult: TColumnMatrix;
  begin

    SetLength(finalResult, configuration.size[dimension xor 1]);

    if dimension = 0 then
      begin
        columnSize := configuration.size[dimension xor 1];

          last := 0;
          for i := (columnSize * column) to ((columnSize * column) + columnSize - 1) do begin
            finalResult[last] := @matrix.cells[i];
            Inc(last);
          end;
      end
    else
      begin

        last := 0;
        step := configuration.size[dimension];

        for i := column to (Length(matrix.cells) - 1) do begin
          if ((i - column) mod step) = 0 then begin
            finalResult[last] := @matrix.cells[i];
            Inc(last);
          end;
        end;
      end;

      GetMatrixColumn := finalResult;

  end;

  // GETS MATRIX ###############################################
  function GetMatrix(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; dimension: Integer): TDimensionMatrix;
    var finalResult: TDimensionMatrix;
    var i: Integer;
  begin
    SetLength(finalResult, configuration.size[dimension], configuration.size[dimension xor 1]);

    for i := 0 to configuration.size[dimension] - 1 do
      finalResult[i] := GetMatrixColumn(matrix, configuration, dimension, i);

    GetMatrix := finalResult;
  end;

  // OUTPUTS DIMENSION MATRIX #######################################

  procedure OutputDimensionMatrix(var configuration: TInitialConfiguration; var matrix: TDimensionMatrixWrapper; var outputDevice: Text);
    var finalResult: TDimensionMatrix;
    var i, j: Integer;
    var character: String;
  begin

    finalResult := GetMatrix(matrix, configuration, 1);

    // Vypisuje
//  for i := 0 to configuration.size[1] - 1 do
//    begin
//      for j := 0 to configuration.size[0] - 1 do
//      begin
//        character := CellMatrixToChr(finalResult[i][j].state);
//        character := ' ' + IntToStr(finalResult[i][j].blocks[1]);
//        if Length(character) = 2 then
//          character := ' ' + character;
//
//        if j = configuration.size[0] - 1 then
//          WriteLn(character)
//        else
//          Write(character);
//      end;
//    end;

    finalResult := GetMatrix(matrix, configuration, 1);

    // Vypisuje
    for i := 0 to configuration.size[1] - 1 do
    begin
      for j := 0 to configuration.size[0] - 1 do
      begin
       character := CellMatrixToChr(finalResult[i][j].state);

        if j = configuration.size[0] - 1 then
          WriteLn(outputDevice, character)
        else
          Write(outputDevice, character);
      end;
    end;

  end;

  // BYTE SET LENGTH ###############################################
  function ByteSetLength(var targetSet: TByteSet; maxLength: Integer): Byte;
    var i: Byte;
    var length: Byte;
  begin
    length := 0;
    for i := 0 to maxLength do
      if i in targetSet then
        Inc(length);
    ByteSetLength := length;
  end;

  // BYTE SET TO ARRAY #############################################
  function ByteSetToArray(var convertedSet: TByteSet; maxLength: Integer): T1DByteArray;
    var i, position: Byte;
    var finalResult: T1DByteArray;
  begin
    SetLength(finalResult, ByteSetLength(convertedSet, maxLength));
    position := 0;
    for i := 0 to maxLength do begin
      if i in convertedSet then begin
        finalResult[position] := i;
        Inc(position);
      end;
    end;
    ByteSetToArray := finalResult;
  end;

  // GET BYTE SET ITEM #############################################
  function GetByteSetItem(var convertedSet: TByteSet; maxLength: Integer; position: Integer): Integer;
    var reachedPosition: Integer;
    var i: Integer;
    var finalResult: Integer;
  begin
    i := 0;
    reachedPosition := -1;

    while (reachedPosition < position) and (i < maxLength) do begin
      if i in convertedSet then begin
        Inc(reachedPosition);
        Break;
      end;
      Inc(i);
    end;
    if reachedPosition < position then
      finalResult := -1
    else
      finalResult := i;

    GetByteSetItem := finalResult;
  end;

  // BYTE SET MIN ################################################
  function ByteSetMin(var targetSet: TByteSet; maxLength: Integer): Integer;
  begin
    ByteSetMin := GetByteSetItem(targetSet, maxLength, 0);
  end;

  // BYTE SET MAX ################################################
  function ByteSetMax(var targetSet: TByteSet; maxLength: Integer): Integer;
    var i: Integer;
    var finalResult: Integer;
  begin
    finalResult := -1;
    for i := maxLength downto 0 do begin
      if i in targetSet then begin
        finalResult := i;
        Break;
      end;
    end;
    ByteSetMax := finalResult;
  end;

  // FINALIZE BLOCK ###########################################################
  procedure FinalizeBlock(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; var cell: TCellMatrixWrapperPointer; var block: TColumnDataBlock; dimension: Integer; column: Integer; blockNum: Integer);
    var currentMatrixColumn: TColumnMatrix;
    var i: Integer;
  begin
    currentMatrixColumn := GetMatrixColumn(matrix, configuration, dimension, column);
    for i := block.result[fromCell] to block.result[toCell] do begin
      if not ((currentMatrixColumn[i].state = cellFinal) or (currentMatrixColumn[i].state = cellFull)) then
        Inc(configuration.final);
      currentMatrixColumn[i].state := cellFinal;
    end;

    // Vyrazuje blok z indexu
    for i := 0 to block.result[fromCell] - 1 do
      Exclude(currentMatrixColumn[i].possibleBlocks[dimension], blockNum);
    for i := block.result[toCell] + 1 to configuration.size[dimension xor 1] - 1 do
      Exclude(currentMatrixColumn[i].possibleBlocks[dimension], blockNum);

  end;

  // UPDATE BLOCK INFO ########################################################
  procedure UpdateBlockInfo(var block: TColumnDataBlock; cellNum: Integer; increase: Boolean);
  begin

    if block.result[fromCell] = -1 then
      begin
        block.result[fromCell] := cellNum;
        block.result[toCell] := cellNum;
      end
    else if cellNum < block.result[fromCell] then
      block.result[fromCell] := cellNum
    else if cellNum > block.result[toCell] then
      block.result[toCell] := cellNum;

    if increase then
      Inc(block.final);

  end;

  // SETS MATRIX CELL BLOCK ###################################################
  // Jako vstupni hodnotu NEOCEKAVA blok -1!
  procedure SetMatrixCellBlock(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; var cell: TCellMatrixWrapperPointer; dimension: Integer; column: Integer; cellNum: Integer; block: Integer);
    var i: Integer;
    var neighbourCell: TCellMatrixWrapperPointer;
    var oldBlock: Integer;
    var blockData: ^TColumnDataBlock;
    var blocks: TBlockList;
  begin

    // Nastavuje blok
    oldBlock := cell.blocks[dimension];
    cell.blocks[dimension] := block;
    cell.possibleBlocks[dimension] := [block];

    blockData := @configuration.matrix[dimension][column].data[block];

    // Nastavuje informace bloku a bunek v konfiguraci
    if not (oldBlock = block) then
      Inc(blockData.final);

    blockData.state := blockStarted;
    UpdateBlockInfo(blockData^, cellNum, False);


    // Nastavuje blok vsem blokum ktere jej nastaveny nemaji doleva i doprava
    //  v ramci teto dimenze
    for i := (cellNum - 1) downto 0 do begin
      neighbourCell := GetMatrixCell(matrix, configuration, dimension, column, i);
      if neighbourCell.state = cellFull then
        begin
         if neighbourCell.blocks[dimension] = -1 then
          begin
            neighbourCell.blocks[dimension] := block;
            UpdateBlockInfo(blockData^, i, True);
          end
         else
            Break;
        end
      else
        begin
          if (cell.state = cellFull) and (neighbourCell.state = cellUnknown) then
            neighbourCell.possibleBlocks[dimension] := neighbourCell.possibleBlocks[dimension] * [block];
          Break;
        end;
    end;
    for i := (cellNum + 1) to (configuration.size[dimension xor 1] - 1) do begin
      neighbourCell := GetMatrixCell(matrix, configuration, dimension, column, i);
      if neighbourCell.state = cellFull then
        begin
          if neighbourCell.blocks[dimension] = -1 then
            begin
              neighbourCell.blocks[dimension] := block;
              UpdateBlockInfo(blockData^, i, True);
            end
          else
            Break;
        end
      else
        begin
          if (cell.state = cellFull) and (neighbourCell.state = cellUnknown) then
            neighbourCell.possibleBlocks[dimension] := neighbourCell.possibleBlocks[dimension] * [block];
          Break;
        end;
    end;


    // Odstavuje vsem nizsim/vyssim bunkam vyssi/nizsi bloky

    if ((cell.blocks[dimension] + 1) <= (blockData.length - 1)) and ((blockData.result[fromCell] - 2) >= 0) then begin
      blocks := [(cell.blocks[dimension] + 1)..(blockData.length - 1)];
      for i := 0 to (blockData.result[fromCell] - 2) do begin
        neighbourCell := GetMatrixCell(matrix, configuration, dimension, column, i);
        neighbourCell.possibleBlocks[dimension] := neighbourCell.possibleBlocks[dimension] - blocks;
      end;
    end;
    if ((cell.blocks[dimension] - 1) >= 0) and ((configuration.size[dimension xor 1] - 1) >= (blockData.result[toCell] + 2)) then begin
      blocks := [0..(cell.blocks[dimension] - 1)];
      for i := (blockData.result[toCell] + 2) to (configuration.size[dimension xor 1] - 1) do begin
        neighbourCell := GetMatrixCell(matrix, configuration, dimension, column, i);
        neighbourCell.possibleBlocks[dimension] := neighbourCell.possibleBlocks[dimension] - blocks;
      end;
    end;

    if blockData.final = blockData.length then
      FinalizeBlock(matrix, configuration, cell, blockData^, dimension, column, block);

  end;

  // SETS MATRIX CELL #########################################################
  procedure SetMatrixCellFull(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; var cell: TCellMatrixWrapperPointer; dimension: Integer; column: Integer; cellNum: Integer; block: Integer);
    var currentConfigColumnLength: Integer;
    var i: Integer;
    var count, matrixLength: Integer;
//    var otherBlocksArray: TColumnData;
//    var firstBlock: Integer;
//    var newBlocks, oldBlocks: TBlockList;
    var otherBlocks: TBlockList;
//    var done: Boolean;
    var workingColumn, workingCell: Integer;
    var neighbourCell: TCellMatrixWrapperPointer;
    var neighbourExists: Boolean;
//    var maxBlock: Integer;
//    var length: Integer;
//    var currentConfigColumn: ^TColumnData;
  begin

    if cell.state = cellUnknown then begin
      cell.state := cellFull;
      Inc(configuration.final);
    end;

//    if state = cellFull then begin

      if not (block = -1) then
        SetMatrixCellBlock(matrix, configuration, cell, dimension, column, cellNum, block);

      for i := 0 to configuration.dimensions - 1 do begin
        if cell.blocks[i] = -1 then
          begin

            if i = (dimension xor 1) then
              begin
                workingColumn := cellNum;
                workingCell := column;
              end
            else
              begin
                workingColumn := column;
                workingCell := cellNum;
              end;

//            currentConfigColumn := @configuration.matrix[i][workingColumn].data;
            currentConfigColumnLength := configuration.matrix[i][workingColumn].length;
            matrixLength := configuration.size[i xor 1];
            count := ByteSetLength(cell.possibleBlocks[i], currentConfigColumnLength);

            if count = 0 then
              begin
//                write('');
                // ERROR, pokus o nastaveni bunky do mista kde se v radku ci sloupci nemohou vyskytnou zadne bloky
                raise Exception.Create('Pokus o nastaveni bunky do mista kde se v radku ci sloupci nemohou vyskytnou zadne bloky. [DIMENSION: ' + IntToStr(dimension) + ', COLUMN: ' + IntToStr(column) + ', CELL: '  + IntToStr(cellNum) + ']. Kontaktujte autora programu')
              end
            else if count = 1 then
              SetMatrixCellBlock(matrix, configuration, cell, i, workingColumn, workingCell, GetByteSetItem(cell.possibleBlocks[i], currentConfigColumnLength, 0))
            else
              begin

                // Dobre, zkusime se rozhlednout doleva a doprava a analyzovat v cem asi tak muzeme byt
                otherBlocks := cell.possibleBlocks[i];
                neighbourExists := False;

                if workingCell - 1 >= 0 then begin
                  neighbourCell := GetMatrixCell(matrix, configuration, i, workingColumn, workingCell - 1);
                  if neighbourCell.state = cellFull then begin
                    otherBlocks := otherBlocks * neighbourCell.possibleBlocks[i];
                    neighbourExists := True;
                  end;
                end;

                if workingCell + 1 < matrixLength then begin
                  neighbourCell := GetMatrixCell(matrix, configuration, i, workingColumn, workingCell + 1);
                  if neighbourCell.state = cellFull then begin
                    otherBlocks := otherBlocks * neighbourCell.possibleBlocks[i];
                    neighbourExists := True;
                  end;
                end;

                count := ByteSetLength(otherBlocks, currentConfigColumnLength);
                if (count = 0) and neighbourExists then
                  begin
                    // ERROR, pokus o nastaveni bunky kde zadna byt nemuze protoze by spojovala dva rozdilne bloky
                    raise Exception.Create('Pokus o nastaveni bunky kde zadna byt nemuze protoze by spojovala dva rozdilne bloky. [DIMENSION: ' + IntToStr(dimension) + ', COLUMN: ' + IntToStr(column) + ', CELL: '  + IntToStr(cellNum) + ']. Kontaktujte autora programu')
                  end
                else if count = 1 then
                  SetMatrixCellBlock(matrix, configuration, cell, i, workingColumn, workingCell, GetByteSetItem(otherBlocks, currentConfigColumnLength, 0))
                else
                  begin end;    // ZATIM TEDY NIC, NO
//                   begin
//
//                     // No nic, z okoli moc moudri nejsme, zkusime se podivat kolik mame mista k zacatku/ke konci a podle toho
//                     // zkusit zjistit nebo alespon eliminovat bloky.
//                     maxBlock := ByteSetMax(cell.possibleBlocks[i], currentConfigColumnLength);
//                     length := 0;
//
//                     for j := 0 to maxBlock do begin
//                       if length > workingCell then
//                         Break;
//                       length := length + currentConfigColumn^[j].length + 1;
//                     end;
//
//                     if  j < maxBlock then begin
//                       oldBlocks := cell.possibleBlocks[i];
//                       newBlocks := [0..j];
//                       cell.possibleBlocks[i] := cell.possibleBlocks[i] * newBlocks;
//
//                       if ByteSetLength(cell.possibleBlocks[i], maxBlock) = 1 then
//                         SetMatrixCellBlock(matrix, configuration, cell, i, column, cellNum, GetByteSetItem(cell.possibleBlocks[i], currentConfigColumnLength, 0));
//
//                       // Sousednim bunkam nazpet dela prunik
//                       for k := workingCell - 1 downto 0 do begin
//                         neighbourCell := GetMatrixCell(matrix, configuration, i, workingColumn, k);
//                         neighbourCell.possibleBlocks[i] := neighbourCell.possibleBlocks[i] * newBlocks;
//                       end;
//
//                     end;
//                  end;
              end;
          end;
//      end;
    end;
  end;

  // SET MATRIX CELL EMPTY ###############################################
  procedure SetMatrixCellEmpty(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; var cell: TCellMatrixWrapper; dimension: Integer; column: Integer; cellNum: Integer);
    var i, j, k, l, m: Integer;
    var blocks: T1DByteArray;
    var workingColumn, workingCell: Integer;
    var secondDimension: Integer;
    var currentMatrixColumn: TColumnMatrix;
    var currentConfigColumn: ^TColumnData;
    var currentBlock: Byte;
    var counter: Integer;
    var forceClose: Byte;
//    blockData := @configuration.matrix[dimension][column].data[block];
//    var blockData: ^TColumnDataBlock;
  begin
    cell.state := cellEmpty;
    secondDimension := dimension xor 1;

    for i := 0 to (configuration.dimensions - 1) do begin
      if i = secondDimension then
        begin
          workingColumn := cellNum;
          workingCell := column;
        end
      else
        begin
          workingColumn := column;
          workingCell := cellNum;
        end;


      currentMatrixColumn := GetMatrixColumn(matrix, configuration, i, workingColumn);
      currentConfigColumn := @configuration.matrix[i][workingColumn].data;

      // Analyzuje jak prazdna bunka rozrizla potencialni bloky a podle toho upravuje index
      if configuration.matrix[i][workingColumn].length > 0 then begin
        blocks := ByteSetToArray(cell.possibleBlocks[i], configuration.matrix[i][workingColumn].length - 1);

        for j := 0 to Length(blocks) - 1 do begin

          currentBlock := blocks[j];
          forceClose := 0;
          k := 0;
          l := 0;

          if workingCell - 1 >= 1 then begin

            counter := 0;
            for k := (workingCell - 1) downto 0 do begin
              if (currentMatrixColumn[k].state = cellEmpty) or not (currentBlock in currentMatrixColumn[k].possibleBlocks[i]) then
                Break
              else
                Inc(counter);

              if (currentMatrixColumn[k].state = cellFull) and (currentMatrixColumn[k].blocks[i] = currentBlock) then
                forceClose := forceClose or 2;
            end;

            Inc(k);

            if (counter > 0) and (currentConfigColumn^[currentBlock].length > counter) then
              forceClose := forceClose or 1;

          end;


          if (workingCell + 1) < configuration.size[i xor 1] then begin

            counter := 0;
            for l := (workingCell + 1) to (configuration.size[i xor 1] - 1) do begin
              if (currentMatrixColumn[l].state = cellEmpty) or not (currentBlock in currentMatrixColumn[l].possibleBlocks[i]) then
                Break
              else
                Inc(counter);

              if (currentMatrixColumn[l].state = cellFull) and (currentMatrixColumn[l].blocks[i] = currentBlock) then
                forceClose := forceClose or 1;

            end;

            Dec(l);

            if (counter > 0) and (currentConfigColumn^[currentBlock].length > counter) then
              forceClose := forceClose or 2;

          end;

          // Uzavira co je k uzavreni
          if (forceClose and 1) = 1 then begin
//            writeln(k);
//            if k = 0 then
//              k := 1;
            for m := k to (workingCell - 1) do begin
              Exclude(currentMatrixColumn[m].possibleBlocks[i], currentBlock);
//              currentMatrixColumn[m].state := cellFinal;
            end;
          end;
          if (forceClose and 2) = 2 then begin
//            if l = configuration.size[i xor 1] - 1 then
//              Dec(l);
            for m := (workingCell + 1) to l do
              Exclude(currentMatrixColumn[m].possibleBlocks[i], currentBlock);
          end;

        end;
      end;

      // Nuluje mozne bloky
      cell.possibleBlocks[i] := [];
    end;
    Inc(configuration.final);
  end;

  // FINDS EMPTY CELLS ###################################################
  procedure FindEmptyCells(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; dimension: Integer; column: Integer);
    var i: Integer;
//    var cell: ^TCellMatrixWrapper;
    var currentMatrixColumn: TColumnMatrix;
  begin
    currentMatrixColumn := GetMatrixColumn(matrix, configuration, dimension, column);
    for i := 0 to Length(currentMatrixColumn) - 1 do
      if (currentMatrixColumn[i].state = cellUnknown) and (currentMatrixColumn[i].possibleBlocks[dimension] = []) then
        SetMatrixCellEmpty(matrix, configuration, currentMatrixColumn[i]^, dimension, column, i);

//    for i := 0 to Length(matrix.cells) - 1 do begin
//      cell := @matrix.cells[i];
//      if cell.state = cellUnknown then begin
//        for j := 0 to configuration.dimensions - 1 do begin
//          if cell.possibleBlocks[j] = [] then begin
//             SetMatrixCellEmpty(cell^, configuration);
//             Break;
//          end;
//        end;
//      end;
//    end;
  end;


  // ANALYZES BLOCK POSITIONS ##################################################
  procedure AnalyzeBlockPositions(var configuration: TInitialConfiguration; var matrix: TDimensionMatrixWrapper; dimension: Integer; column: Integer);
    var i, j: Integer;
    var secondDimension: Integer;
    var matrixLength: Integer;
    var currentMatrixColumn: TColumnMatrix;
    var currentConfigColumn: ^TColumnData;
    var leftBorder, rightBorder: Integer;
  begin

    secondDimension := dimension xor 1;
    matrixLength := configuration.size[secondDimension];
    currentConfigColumn := @configuration.matrix[dimension][column].data;
    currentMatrixColumn := GetMatrixColumn(matrix, configuration, dimension, column);

    // ANALYZA VSECH MOZNYCH POZIC VSECH BLOKU

    // Pro vsechny bloky vetsi nez delta mezi vynucenymi bunkami a celkovou delkou sloupce
    //   analyzuje jimi v kazdem pripade obsazene bunky.
    for i := 0 to configuration.matrix[dimension][column].length - 1 do begin
      if not (currentConfigColumn^[i].state = blockFinished) then begin

          // from left
          for leftBorder := 0 to (matrixLength - 1) do
            if not (currentMatrixColumn[leftBorder].state = cellEmpty) then
              Break;

          for j := 0 to i - 1 do begin
            leftBorder := leftBorder + currentConfigColumn^[j].length;
            Inc(leftBorder);
          end;


          for rightBorder := (matrixLength - 1) downto 0 do
            if not (currentMatrixColumn[rightBorder].state = cellEmpty) then
              Break;

          for j := i + 1 to matrixLength do begin
            if currentConfigColumn^[j].length = 0 then
              Break
            else
              begin
                rightBorder := rightBorder - currentConfigColumn^[j].length;
                Dec(rightBorder);
              end
          end;

          // Index
//          if dimension = 1 then begin
//            WriteLn(IntToStr(dimension) + ':' + IntToStr(column));
//            writeln(IntToStr(i) + ' ' + IntToStr(leftBorder) + ' to ' + IntToStr(rightBorder));
//          end;

          for j := leftBorder to rightBorder do
            Include(currentMatrixColumn[j].possibleBlocks[dimension], i);

          currentConfigColumn^[i].range[fromCell] := leftBorder;
          currentConfigColumn^[i].range[toCell] := rightBorder;

        end;
      end;
  end;

  // REANALYZES BLOCK POSITIONS ################################################
  procedure ReanalyzeBlockPositions(var configuration: TInitialConfiguration; var matrix: TDimensionMatrixWrapper; dimension: Integer; column: Integer);
    var i, j: Integer;
//    var secondDimension: Integer;
//    var matrixLength: Integer;
    var currentMatrixColumn: TColumnMatrix;
    var currentConfigColumn: ^TColumnData;
//    var leftBorder, rightBorder: Integer;
    var skip: Boolean;
    var currentConfigBlock: ^TColumnDataBlock;
  begin
//    secondDimension := dimension xor 1;
//    matrixLength := configuration.size[secondDimension];
    currentConfigColumn := @configuration.matrix[dimension][column].data;
    currentMatrixColumn := GetMatrixColumn(matrix, configuration, dimension, column);


    // ANALYZA VSECH MOZNYCH POZIC VSECH BLOKU

    // Pro vsechny bloky vetsi nez delta mezi vynucenymi bunkami a celkovou delkou sloupce
    //   analyzuje jimi v kazdem pripade obsazene bunky.
//    if (dimension = 0) and (column = 23) then
//      write('');
    for i := 0 to configuration.matrix[dimension][column].length - 1 do begin
      currentConfigBlock := @currentConfigColumn^[i]; 
      if not (currentConfigBlock.state = blockFinished) and (currentConfigBlock.range[fromCell] >= 0) and (currentConfigBlock.range[toCell] >= 0)then begin

        for j := currentConfigBlock.range[fromCell] to currentConfigBlock.range[toCell] do
//          if ((currentMatrixColumn[j].state = cellUnknown) and (i in currentMatrixColumn[j].possibleBlocks[dimension])) or (((currentMatrixColumn[j].state = cellFull) or (currentMatrixColumn[j].state = cellFinal)) and (currentMatrixColumn[j].blocks[dimension] = i)) then
          if i in currentMatrixColumn[j].possibleBlocks[dimension] then
            Break;
        if currentConfigBlock.range[toCell] - (j - 1) < currentConfigBlock.length then
          // ERROR, mozny vyskyt bloku se zkratil natolik, ze se uz by se sam do sebe nevesel
          raise Exception.Create('Mozny vyskyt bloku se zkratil natolik, ze se uz by se sam do sebe nevesel. [DIMENSION: ' + IntToStr(dimension) + ', COLUMN: ' + IntToStr(column) + ', BLOCK: '  + IntToStr(i) + ']. Kontaktujte autora programu')
        else
          currentConfigBlock.range[fromCell] := j;
//        if currentConfigColumn^[i].range[fromCell] > currentConfigColumn^[i].range[toCell] then
//          write('');
        for j := currentConfigBlock.range[toCell] downto currentConfigBlock.range[fromCell] do
//          if ((currentMatrixColumn[j].state = cellUnknown) and (i in currentMatrixColumn[j].possibleBlocks[dimension])) or (((currentMatrixColumn[j].state = cellFull) or (currentMatrixColumn[j].state = cellFinal)) and (currentMatrixColumn[j].blocks[dimension] = i)) then
          if i in currentMatrixColumn[j].possibleBlocks[dimension] then
            Break;
        if currentConfigBlock.range[fromCell] + (j + 1) < currentConfigBlock.length then
          // ERROR, mozny vyskyt bloku se zkratil natolik, ze se uz by se sam do sebe nevesel 
          raise Exception.Create('Mozny vyskyt bloku se zkratil natolik, ze se uz by se sam do sebe nevesel. [DIMENSION: ' + IntToStr(dimension) + ', COLUMN: ' + IntToStr(column) + ', BLOCK: '  + IntToStr(i) + ']. Kontaktujte autora programu')
        else
          currentConfigBlock.range[toCell] := j;

      end;
    end;

    // JSOU-LI VSECHNY BLOKY V UVODNIM STAVU, PREPOCITA HRANICE PODLE KRAJNICH PRAZDNYCH BUNEK

    skip := False;    
    for i := 0 to configuration.matrix[dimension][column].length - 1 do
      if not (currentConfigColumn^[i].state = blockInitial) then
        skip := True;
    for i := currentConfigColumn^[0].range[fromCell] to currentConfigColumn^[configuration.matrix[dimension][column].length - 1].range[toCell] do
//      if (currentMatrixColumn[i].state = cellEmpty) or ((currentMatrixColumn[i].state = cellFinal) or (currentMatrixColumn[i].state = cellFull) and (currentMatrixColumn[i].blocks[dimension] <= 0)) then
      if currentMatrixColumn[i].state = cellEmpty then
        skip := True;      

    if skip = False then
      AnalyzeBlockPositions(configuration, matrix, dimension, column);

  end;

  // COMPUTES COLUMN - METHOD 1 ###############################################

  procedure ComputeColumn_Method1(dimension: Integer; column: Integer; var configuration: TInitialConfiguration; var matrix: TDimensionMatrixWrapper);
    var i, j: Integer;
    var currentMatrixColumn: TColumnMatrix;
    var currentConfigColumn: TColumnData;
    var currentColumn: ^TColumnDataBlock;
    var fromLeft, fromRight: Integer;
  begin

    currentConfigColumn := configuration.matrix[dimension][column].data;
    currentMatrixColumn := GetMatrixColumn(matrix, configuration, dimension, column);
//          if (dimension = 0) and (column = 23) then
//            write('');

    // ANALYZA RADKU/SLOUPCE NASTAVENIM VSECH NEZAHAJENYCH BLOKU DO EXTREMNICH POZIC
    for i := 0 to configuration.matrix[dimension][column].length - 1 do begin

      currentColumn := @currentConfigColumn[i];

      if currentColumn.state = blockInitial then begin

        fromRight := currentColumn.range[toCell] - currentColumn.length + 1;
        fromLeft := currentColumn.range[fromCell] + currentColumn.length - 1;

        // Intersection
        for j := fromRight to fromLeft do begin
          SetMatrixCellFull(matrix, configuration, currentMatrixColumn[j], dimension, column, j, i);
//          WriteLn(IntToStr(dimension) + ':' + IntToStr(column) + ':' + IntToStr(j));
//          Exit;
        end;

      end;
    end;
  end;

  procedure ComputeColumn_Method2(dimension: Integer; column: Integer; var configuration: TInitialConfiguration; var matrix: TDimensionMatrixWrapper);
    var i, j: Integer;
    var secondDimension: Integer;
    var matrixLength: Integer;
    var currentMatrixColumn: TColumnMatrix;
    var currentConfigColumn: TColumnData;
    var border, originalBorder: Integer;
    var fromLeft, fromRight: Integer;
    var blockLength, currentBlockLength: Integer;
    var currentBlock: Integer;
    var skip: Byte;
  begin

//    WriteLn(IntToStr(dimension) + ':' + IntToStr(column));

    secondDimension := dimension xor 1;
    matrixLength := configuration.size[secondDimension];
    currentConfigColumn := configuration.matrix[dimension][column].data;
    currentMatrixColumn := GetMatrixColumn(matrix, configuration, dimension, column);

    // ANALYZA RADKU/SLOUPCE POSOUVANIM BLOKU VUCI Z NEJ ZNAMYM POLICKUM DO EXTREMNICH POZIC

    i := 0;
    while i < matrixLength do begin
      if (currentMatrixColumn[i].state = cellFull) or (currentMatrixColumn[i].state = cellFinal) then begin
        currentBlock := currentMatrixColumn[i].blocks[dimension];
        if (currentBlock > -1) and (currentConfigColumn[currentBlock].state = blockStarted) then begin
          blockLength := 0;
          currentBlockLength := currentConfigColumn[currentBlock].length;
          while ((i + blockLength) < matrixLength) and ((currentMatrixColumn[i + blockLength].state = cellFull) or (currentMatrixColumn[i + blockLength].state = cellFinal)) do
            Inc(blockLength);                              // Pocita delku sekvence plnych bunek


          skip := 0;

          // Vyhledava pravou pozici umisteni co nejvice doleva
//          if (dimension =0) and (column = 6) then
//            write('');
          border := i + blockLength - currentBlockLength;
          if border < 0 then
            begin
              fromLeft := 0;
              border := 0;
            end
          else
            fromLeft := border;
          originalBorder := border;
//          for k := border to (i - 1) do
//            if not (currentBlock in currentMatrixColumn[k].possibleBlocks[dimension]) then
//              border := k;
          while (border <= (matrixLength - 1)) and (not (currentBlock in currentMatrixColumn[border].possibleBlocks[dimension])) and (not (border = i)) do
            Inc(border);

          if border >= matrixLength then
            skip := skip or 1;

          // Vice doleva odstranuje tento blok z indexu bloku
          if (skip and 1) = 0 then begin
            fromLeft := fromLeft + (border - originalBorder);
            for j := 0 to (fromLeft - 1) do
              Exclude(currentMatrixColumn[j].possibleBlocks[dimension], currentBlock);
            currentConfigColumn[currentBlock].range[fromCell] := (fromLeft);// - 1);
//            if fromleft = 0 then
//              write('');

            fromLeft := fromLeft + currentBlockLength - 1;
          end;

          // Vyhledava levou pozici umisteni co nejvice doprava
//          if (dimension = 1) and (column = 8) and (i = 10) then
//            write('');
          border := i + currentBlockLength - 1;
          if border > (matrixLength - 1) then
            begin
              fromRight := (matrixLength - 1);
              border := (matrixLength - 1);
            end
          else
            fromRight := border;
          originalBorder := border;
//          for k := border downto (i + 1) do
//            if not (currentBlock in currentMatrixColumn[k].possibleBlocks[dimension]) then
//              border := k;
          while (border >= 0) and (not (currentBlock in currentMatrixColumn[border].possibleBlocks[dimension])) and (not (border = i)) do
            Dec(border);

          if border = 0 then begin
//            border := 0;
            skip := skip or 2;
          end;

          // Vice doprava odstranuje tento blok z indexu bloku
          if (skip and 2) = 0 then begin
            fromRight := fromRight - (originalBorder - border);
            for j := (fromRight + 1) to (matrixLength - 1) do
              Exclude(currentMatrixColumn[j].possibleBlocks[dimension], currentBlock);
//            if currentConfigColumn[currentBlock].range[fromCell] > fromRight then
//              write('');
            currentConfigColumn[currentBlock].range[toCell] := (fromRight);// + 1);
//            if currentConfigColumn[currentBlock].range[fromCell] > currentConfigColumn[currentBlock].range[toCell] then
//              write('');

            fromRight := fromRight - currentBlockLength + 1;
          end;

          // ####

          if skip = 0 then
            begin
              // Prunik
//              if (dimension = 1) and (column = 8) and (i = 10) then
//                write('');
              if (fromRight < 0) or (fromLeft >= matrixLength) then 
                // ERROR, pokus o nastaveni neexistujici bunky
                raise Exception.Create('Pokus o nastaveni neexistujici bunky. [METODA: 2, DIMENSION: ' + IntToStr(dimension) + ', COLUMN: ' + IntToStr(column) + ', FROM LEFT: '  + IntToStr(fromLeft) + ', FROM RIGHT: '  + IntToStr(fromRight) + ']. Kontaktujte autora programu');
              for j := fromRight to fromLeft do begin
//                writeln(IntToStr(dimension) + ':' + IntToStr(column) + ':' + IntToStr(i));
                  SetMatrixCellFull(matrix, configuration, currentMatrixColumn[j], dimension, column, j, currentBlock);
              end;

              if fromRight < fromLeft then
                i := fromLeft + 1
              else
                i := i + blockLength;

//              WriteLn(IntToStr(currentBlock) + ' ' + IntToStr(fromleft) + ' to ' + IntToStr(fromRight));
            end
          else
            begin
              Inc(i);
              Continue;
            end;

        end;

      end;
      Inc(i);
    end;

  end;

  // GLUES PIECES ##################################################

  procedure GluePieces(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; dimension: Integer; column: Integer);
    var i, j: Integer;
    var currentMatrixColumn: TColumnMatrix;
    var currentConfigColumn: TColumnDataWrapper;
    var currentConfigBlock: TColumnDataBlock;
  begin
    currentConfigColumn := configuration.matrix[dimension][column];
    currentMatrixColumn := GetMatrixColumn(matrix, configuration, dimension, column);

    // Pro vsechny bloky
    for i := 0 to currentConfigColumn.length do begin
      currentConfigBlock := currentConfigColumn.data[i];
      if (currentConfigBlock.result[fromCell] >= 0) and (currentConfigBlock.result[toCell] > currentConfigBlock.result[fromCell]) then
        for j := currentConfigBlock.result[fromCell] to currentConfigBlock.result[toCell] do
          if not ((currentMatrixColumn[j].state = cellFinal) or (currentMatrixColumn[j].state = cellFull)) then
            SetMatrixCellFull(matrix, configuration, currentMatrixColumn[j], dimension, column, j, i);
    end;
  end;

  // CHECKS FOR FINAL CELL POSITIONS ###############################
  procedure CheckCellsForFinal(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration; dimension: Integer; column: Integer);
    var currentMatrixColumn: TColumnMatrix;
    var secondDimension: Integer;
    var i, j: Integer;
    var currentConfigColumnLength: Integer;
    var cell: TCellMatrixWrapperPointer;
    var workingColumn, workingCell: Integer;
    var count: Integer;
  begin
    secondDimension := dimension xor 1;
    currentMatrixColumn := GetMatrixColumn(matrix, configuration, dimension, column);
    for i := 0 to (configuration.size[secondDimension] - 1) do begin
      cell := currentMatrixColumn[i];
      if cell.state = cellFull then begin

        for j := 0 to (configuration.dimensions - 1) do begin
        if cell.blocks[j] = -1 then begin

            if j = secondDimension then
              begin
                workingColumn := i;
                workingCell := column;
              end
            else
              begin
                workingColumn := column;
                workingCell := i;
              end;

            currentConfigColumnLength := configuration.matrix[j][workingColumn].length;
            count := ByteSetLength(cell.possibleBlocks[j], currentConfigColumnLength);

            if count = 1 then
              SetMatrixCellBlock(matrix, configuration, cell, j, workingColumn, workingCell, GetByteSetItem(cell.possibleBlocks[j], currentConfigColumnLength, 0));

          end;
        end;
      end;
    end;
  end;

  // INITIALIZE ####################################################

  procedure Initialize(var configuration: TInitialConfiguration; var matrix: TDimensionMatrixWrapper);
    var i, j: Integer;
    var cellCount, dimensionCount: Integer;
  begin

    cellCount := 1;
    dimensionCount := configuration.dimensions;

    for i := 0 to dimensionCount - 1 do
      cellCount := cellCount * configuration.size[i];

    SetLength(matrix.cells, cellCount);

    for i := 0 to cellCount - 1 do begin
      matrix.cells[i].state := cellUnknown;
      SetLength(matrix.cells[i].possibleBlocks, dimensionCount);
      SetLength(matrix.cells[i].blocks, dimensionCount);
      for j := 0 to dimensionCount - 1 do begin
        matrix.cells[i].blocks[j] := -1;
        matrix.cells[i].possibleBlocks[j] := [];
      end;
    end;

  end;

  // ITERATES ######################################################
  procedure Iteration(var matrix: TDimensionMatrixWrapper; var configuration: TInitialConfiguration);
    var i, j: Integer;
  begin
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        ComputeColumn_Method1(i, j, configuration, matrix);
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        GluePieces(matrix, configuration, i, j);
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        CheckCellsForFinal(matrix, configuration, i, j);        
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        FindEmptyCells(matrix, configuration, i, j);
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        ReanalyzeBlockPositions(configuration, matrix, i, j);
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        ComputeColumn_Method2(i, j, configuration, matrix);
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        GluePieces(matrix, configuration, i, j);
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        CheckCellsForFinal(matrix, configuration, i, j);
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        FindEmptyCells(matrix, configuration, i, j);
    for i := 0 to configuration.dimensions - 1 do
      for j := 0 to configuration.size[i] - 1 do
        ReanalyzeBlockPositions(configuration, matrix, i, j);
  end;

  // VERSION INFO ##################################################  
  procedure VersionInfo();
  begin
    WriteLn('');  
    WriteLn('Nonogram Solver');
    WriteLn('alpha version, 20080610');
    WriteLn('(c) 2008 Martin Kozak (martinkozak@martinkozak.net)');
    WriteLn('');    
    WriteLn('All rights granted to:');
    WriteLn('Faculty of Nuclear Sciences and Physical Engineering');
    WriteLn('Czech Technical University in Prague');
    WriteLn('');
    WriteLn('Distributed under the terms of MIT License.');    
  end;

  // HELP ##########################################################
  procedure Help();
  begin
    WriteLn('Usage:');
    WriteLn('nonogram-solver file [options]...');
    WriteLn('Options:');
    WriteLn('  --help             Display this information');
    WriteLn('  --version          Display version information');
    WriteLn('  --out <file>       Result to specified file');
    WriteLn('  --verbose [<1|2>]  Level of verbosity (nonogram, additional info)');
    WriteLn('  --interactive      Iterate interactivelly');
  end;

  // ###############################################################

  ////////// VARIABLES
  var configuration: TInitialConfiguration;
  var matrix: TDimensionMatrixWrapper;
  var i, j, counter: Integer;
  var last: Array[0..1] of Integer;
  var parameter: String;
  var parameters: TApplicationParameters;
  var parametersCount: Integer;
  var outFile: Text;

begin

  parameters.verbose := 0;
  parameters.input := '';
  parameters.output := '';
  parameters.interactive := False;

  // Loading parameters

  parametersCount := ParamCount;
  i := 0;

  while i <= parametersCount do
    begin
      parameter := ParamStr(i);
      if parameter = '--version' then
        begin
          VersionInfo();
          Exit;
        end
      else if (parameter = '--help') or (parameter = '-h') then
        begin
          Help();
          Exit;
        end
      else if parameter = '--interactive' then
        begin
          parameters.interactive := True;
          parameters.verbose := 2;
        end
      else if parameter = '--verbose' then
        begin
          Inc(i);
          parameter := ParamStr(i);
          if not (Copy(parameter, 0, 2) = '--') and not (parameter = '') then
            parameters.verbose := StrToInt(ParamStr(i))
          else
            parameters.verbose := 1;
        end
      else if parameter = '--out' then
        begin
          Inc(i);
          parameter := ParamStr(i);
          if (Copy(parameter, 0, 2) = '--') or (parameter = '') then
            begin
              WriteLn('WARNING: Out parameter omitted.');
              Dec(i);
            end
          else
            parameters.output := parameter;
        end
      else if i = 1 then
        parameters.input := parameter;
        
      Inc(i);
    end;

  if parameters.input = '' then
    begin
      WriteLn('No input file specified.');
      Halt(1)
    end
  else if not FileExists(parameters.input) then
    begin
      WriteLn('Input file not exists.');
      Halt(2);
    end;


  // Initialization
  configuration := LoadConfiguration(parameters.input);
  Initialize(configuration, matrix);

  // Initial
  for i := 0 to configuration.dimensions - 1 do
    for j := 0 to configuration.size[i] - 1 do
      AnalyzeBlockPositions(configuration, matrix, i, j);

  // Other
  counter := 0;
  last[0] := -1;
  last[1] := -1;

  while (not ((last[0] = last[1]) and (configuration.final = last[0]))) and (not (configuration.final = Length(matrix.cells))) do begin
    Inc(counter);
    Iteration(matrix, configuration);

    if parameters.interactive = True then begin
      if parameters.verbose > 0 then begin
        OutputDimensionMatrix(configuration, matrix, Output);
        WriteLn('');
      end;
      if parameters.verbose > 1 then begin
        WriteLn('iteration #' + IntToStr(counter));
        WriteLn('Finished: ' + IntToStr(Floor((configuration.final / Length(matrix.cells)) * 100)) + '% (' + IntToStr(configuration.final) + '/' + IntToStr(Length(matrix.cells)) + ' total)');
      end;
      ReadLn;
    end;

    last[counter mod 2] := configuration.final;
  end;

  if parameters.interactive = False then begin
    if parameters.verbose > 0 then begin
      for i := 0 to Length(matrix.cells) - 1 do
        if matrix.cells[i].state = cellFull then
           matrix.cells[i].state := cellFinal;
      OutputDimensionMatrix(configuration, matrix, Output);
      WriteLn('');
    end;
    if parameters.verbose > 1 then begin
      WriteLn('iteration #' + IntToStr(counter));
      WriteLn('Finished: ' + IntToStr(Floor((configuration.final / Length(matrix.cells)) * 100)) + '% (' + IntToStr(configuration.final) + '/' + IntToStr(Length(matrix.cells)) + ' total)');      
    end;
  end;

  if not (parameters.output = '') then begin
    Assign(outFile, parameters.output);
    Rewrite(outFile);
    OutputDimensionMatrix(configuration, matrix, outFile);
    Close(outFile);
  end;

end.
