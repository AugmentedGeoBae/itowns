/*
*Quantized Mesh is Triangulated Irregular Network (TIN) mesh format
*
*https://github.com/CesiumGS/quantized-mesh
*
*/



import * as THREE from 'three.js';



const highwaterDecode = (indices) => {

    var arr = [];

    var highest = 0;
    for (var i = 0; i < indices.length; ++i) {
        var code = indices[i];
        arr.push(highest - code);
        if (code === 0) {
            ++highest;
        }
    }
    return arr;
}

const getEdges = (buff, data, startPos, vertices) => {

  edges = {};



  let indicesDecoded = nOctDecode(edgeIndices);

  edges.west
  edges.south
  edges.east
  edges.north

  return edges;

}

function signNotZero (vec) {
  return new THREE.Vector2( vec.x >= 0 ? 1 : -1, vec.y >= 0 ? 1 : -1 )
}

function nOctDecode (encodedVec) {
  let decodedVec = encodedVec.divideScalar(255).multiplyScalar(2).subScalar(1);

  decodedVec = new THREE.Vector3( decodedVec.x, decodedVec.y, 1 - Math.abs(decodedVec.x) - Math.abs(decodedVec.y) )

  if (decodedVec.z < 0) {
    const xy = new THREE.Vector2(decodedVector.x, decodedVector.y)
    const xyAbs = xy.distanceTo(new THREE.Vector2(0, 0))
    const xySign = signNotZero(xy)
    const decodedXy = xySign.multiplyScalar(1 - xyAbs)

    decodedVec.set(decodedXy.x, decodedXy.y, decodedVec.z)
  }

  return decodedVec.normalize()
}

const getNormals = (buff, data, startPos, normalsLen) => {

  //normals are x,y and 2 byte encoded
  let normalArr = new Uint8Array(buff, starPos, normalsLen);

  const view = new DataView(vertexNormalsBuffer)
  const elementsPerEncodedNormal = 2
  const elementsPerNormal = 3
  const vertexNormalsAttributeArray = new Float32Array(vertexData.length)

  for (let position = 0, i = 0; position < vertexNormalsBuffer.byteLength; position += Uint8Array.BYTES_PER_ELEMENT * elementsPerEncodedNormal, i++) {
    const decodedNormal = decodeOct(new THREE.Vector2(
      view.getUint8(position, true),
      view.getUint8(position + Uint8Array.BYTES_PER_ELEMENT, true)
    ))

    vertexNormalsAttributeArray[i * elementsPerNormal] = decodedNormal.x
    vertexNormalsAttributeArray[i * elementsPerNormal + 1] = decodedNormal.y
    vertexNormalsAttributeArray[i * elementsPerNormal + 2] = decodedNormal.z
  }

  //return new THREE.BufferAttribute(vertexNormalsAttributeArray, elementsPerNormal, true)

  let vNormals = vertexNormalsAttributeArray;

  return vNormals;
}

const getWaterMask = (buff, data, startPos, maskLen) => {

  let waterMask = buff.slice(startPos, startPos + maskLen);

  let type;

  if (waterMask.length() > 1) {
    type = "mixed";
  } else {
    if (data.getUint8(waterMask[0]) === 0) {
      type = "land";
    } else {
      type = "water";
    }
  }

  return waterMask, type;
}

const zigZagDecode = (value) => {
    return (value >> 1) ^ (-(value & 1));
}

const getUint32Array = (data, startPos, count) => {
    return new Uint32Array(data.slice(startPos, startPos + 4*count));
}

const getUint16Array = (data, startPos, count) => {
    return new Uint16Array(data.slice(startPos, startPos + 2*count));
}

const getFloat32Array = (data, startPos, count) => {
    return new Float32Array(data.slice(startPos, startPos + 4*count));
}

const getHeader = (data, byteOffset) => {

  return {
      bytes: data.byteLength,
      centerX : data.getFloat64(byteOffset, true), //getFloat64(data, byteOffset),
      centerY : data.getFloat64(byteOffset + 8, true),
      centerZ : data.getFloat64(byteOffset + 16, true),
      minimumHeight : data.getFloat32(byteOffset + 24, true),
      maximumHeight : data.getFloat32(byteOffset + 28, true),
      boundingSphereCenterX : data.getFloat64(byteOffset + 32, true),
      boundingSphereCenterY : data.getFloat64(byteOffset + 40, true),
      boundingSphereCenterZ : data.getFloat64(byteOffset + 48, true),
      boundingSphereRadius : data.getFloat64(byteOffset + 56, true),
      horizonOcclusionPointX : data.getFloat64(byteOffset + 64, true),
      horizonOcclusionPointY : data.getFloat64(byteOffset + 72, true),
      horizonOcclusionPointZ : data.getFloat64(byteOffset + 80, true)
  }

}

function getVertices(uArray, vArray, heightArray, indexArray) {
    var vertices = [];
    for (var i = 0; i < uArray.length; i++) {
        //vertices.push(new THREE.Vector3(uArray[i]/100, vArray[i]/100, heightArray[i]/200));
        vertices.push(uArray[i]/100);
        vertices.push(vArray[i]/100);
        vertices.push(heightArray[i]/200);
    }
    return vertices;
};

function getFaces(uArray, vArray, heightArray, indexArray) {
    var faces = [];
    for (var i = 0; i < indexArray.length; i+=3) {
        faces.push(new THREE.Face3(indexArray[i], indexArray[i+1], indexArray[i+2]));
    }
    return faces;
};

function getGeometry(indexArray, vertices) {
  let geometry = new THREE.BufferGeometry();
  geometry.setIndex(indexArray);
  geometry.setAttribute( 'position', new THREE.Float32BufferAttribute( vertices, 3 ) );

  return geometry;
};

function getExtensions(buff, data, startPos) {

  let exts = [];
  //let byteOffset1 = byteOffset;

  while (startPos < data.byteLength) {

    var extId = data.getUint8(startPos, true); //1 byte (8 bits)
    startPos += Uint8Array.BYTES_PER_ELEMENT;
    var extLength = data.getUint32(startPos, true);
    startPos += Uint32Array.BYTES_PER_ELEMENT;

    startPos += extLength;

    if (extId == 1) {

      exts.push( getNormals(buff, data, startPos, extLength) );

    } else if (extId == 2) {

      exts.push( getWaterMask(buff, data, startPos, extLength) );

    } else if (extId == 4) {
      exts.push( getMetadata(buff, data, startPos, extLength) );
    } else {
      console.warn(`ID ${extId} not recognised, only 1[= vertex normals] and 2[= water mask] are accepted`);
    }

    //getExtension();

  }

  return exts, startPos;

}

export default {
  /** Parse quantized mesh buffer and extract THREE.Scene and batch table
   * @param {ArrayBuffer} buffer - the quantized mesh buffer.
   * @param {Object} options - additional properties.
   * @param {string=} [options.quantizedMeshUpAxis='Y'] - mesh up axis.
   * @param {string} options.urlBase - the base url of the b3dm file (used to fetch textures for the embedded glTF model).
   * @param {boolean=} [options.doNotPatchMaterial='false'] - disable patching material with logarithmic depth buffer support.
   * @param {float} [options.opacity=1.0] - the b3dm opacity.
   * @param {boolean|Material=} [options.overrideMaterials='false'] - override b3dm's embedded glTF materials. If overrideMaterials is a three.js material, it will be the material used to override.
   * @return {Promise} - a promise that resolves with an object containig a THREE.Scene (gltf) and a batch table (batchTable).
   *
   */


  parse(buffer, options) {

    let view = new DataView(buffer);

    let quantizedMeshComponents = {};

    let byteOffset = 0;

    let header = getHeader(view, byteOffset);
    byteOffset += 88;

    quantizedMeshComponents.header = header;

    var vertexCount = view.getUint32(byteOffset, true);
    byteOffset += Uint32Array.BYTES_PER_ELEMENT;

    var uArray = getUint16Array(buffer, byteOffset, vertexCount);
    byteOffset += vertexCount * Uint16Array.BYTES_PER_ELEMENT;

    var vArray = getUint16Array(buffer, byteOffset, vertexCount);
    byteOffset += vertexCount * Uint16Array.BYTES_PER_ELEMENT;

    var heightArray = getUint16Array(buffer, byteOffset, vertexCount);
    byteOffset += vertexCount * Uint16Array.BYTES_PER_ELEMENT;

    let i;
    let u = 0;
    let v = 0;
    let height = 0;

    for (i = 0; i < uArray.length; ++i) {
        u += zigZagDecode(uArray[i]);
        v += zigZagDecode(vArray[i]);
        height += zigZagDecode(heightArray[i]);

        uArray[i] = u;
        vArray[i] = v;
        heightArray[i] = height;
    }

    if (byteOffset % 2 !== 0) {
        byteOffset += (2 - (byteOffset % 2));
    }


    var triangleCount = view.getUint32(byteOffset, true);
    byteOffset += Uint32Array.BYTES_PER_ELEMENT;

    //Immediately following the vertex data is the index data.
    //Indices specify how the vertices are linked together into triangles.
    //If tile has more than 65536 vertices, the tile uses the IndexData32 structure to encode indices.
    //Otherwise, it uses the IndexData16 structure.

    //To enforce proper byte alignment,
    //padding is added before the IndexData to ensure 2 byte alignment for IndexData16 and 4 byte alignment for IndexData32.

    var indices;

    if (vertexCount > 65536) {

      indices = getUint32Array(buffer, byteOffset, triangleCount * 3);
      byteOffset += triangleCount * 3 * 4;

    } else {

      indices = getUint16Array(buffer, byteOffset, triangleCount * 3);
      byteOffset += triangleCount * 3 * 2;

    }

    let indexArray = highwaterDecode(indices);

/////////////////////////////////

    let vertices = getVertices(uArray, vArray, heightArray, indexArray);

    //why does it calculate faces?
    let faces = getFaces(uArray, vArray, heightArray, indexArray);

    quantizedMeshComponents.geometry = getGeometry(indexArray, vertices);

/////////////////////////////////

    // triangleCount * 3 = vertex count?
    quantizedMeshComponents.edges = getEdges(buffer, view, byteOffset, triangleCount * 3);
    //for (var i = 0, i < edgeIndicesArray, i++) {
    //  byteOffest += edgeIndicesArray[i].numBytes;
    //}
    //byteOffset += edgeIndicesArray.west.length +  length

    quantizedMeshComponents.extensions, extEndPos = getExtensions();

    byteOffset += extEndPos;

    if (byteOffset != view.byteLength) {

      console.warn("byteOffset != view.byteLength");

    }

    //return new ThreeQuantizedMeshTile(header, uArray, vArray, heightArray, indexArray);

    //return Promise.resolve(quantizedMeshTile_bufferGeometry);

    console.log(quantizedMeshComponents);

    return Promise.resolve(quantizedMeshComponents);

    //It should return a buffer geometry as part of an object containing vertexes, indices, edges, extensions

  }
};