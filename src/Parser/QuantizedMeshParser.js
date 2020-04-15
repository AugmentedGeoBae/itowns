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

const edgeIndicesDecode = (data, startPos, vertexCount) => {

  edges = {};

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

const getNormals = (data, startPos, normalsLen) => {

  //normals are x,y and 2 byte encoded
  let normalArr = new Uint8Array(data, starPos, normalsLen);

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

const getWaterMask = (data, startPos, maskLen) => {

  let waterMask = data.slice(startPos, startPos + maskLen);

  let type;

  if (waterMask.length() > 1) {
    type = "mixed";
  } else {
    if (getUint8(waterMask[0]) === 0) {
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
const getUint32 = (data, startPos) => {
    return getUint32Array(data, startPos, 1)[0];
}
const getUint16Array = (data, startPos, count) => {
    return new Uint16Array(data.slice(startPos, startPos + 2*count));
}
const getUint16 = (data, startPos) => {
    return getUint16Array(data, startPos, 1)[0];
}
const getFloat64Array = (data, startPos, count) => {
    return new Float64Array(data.slice(startPos, startPos + 8 * count));
}
const getFloat64 = (data, startPos) => {
    return getFloat64Array(data, startPos, 1)[0];
}
const getFloat32Array = (data, startPos, count) => {
    return new Float32Array(data.slice(startPos, startPos + 4*count));
}
const getFloat32 = (data, startPos) => {
    return getFloat32Array(data, startPos, 1)[0];
}

const getHeader = (data, byteOffset) => {

  return {
      bytes: data.byteLength,
      centerX : getFloat64(data, byteOffset),
      centerY : getFloat64(data, byteOffset + 8),
      centerZ : getFloat64(data, byteOffset + 16),
      minimumHeight : getFloat32(data, byteOffset + 24),
      maximumHeight : getFloat32(data, byteOffset + 28),
      boundingSphereCenterX : getFloat64(data, byteOffset + 32),
      boundingSphereCenterY : getFloat64(data, byteOffset + 40),
      boundingSphereCenterZ : getFloat64(data, byteOffset + 48),
      boundingSphereRadius : getFloat64(data, byteOffset + 56),
      horizonOcclusionPointX : getFloat64(data, byteOffset + 64),
      horizonOcclusionPointY : getFloat64(data, byteOffset + 72),
      horizonOcclusionPointZ : getFloat64(data, byteOffset + 80)
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

function getGeometry(index, vertices) {
  let geometry = new THREE.BufferGeometry();
  geometry.setIndex(index);
  geometry.setAttribute( 'position', new THREE.Float32BufferAttribute( vertices, 3 ) );

  return geometry;
};

function getExtension(data, startPos) {

  let exts = [];
  let byteOffset1 = byteOffset;

  while (byteOffset1 < buffer.byteLength) {

    var extId = getUint8(buffer, byteOffset); //1 byte (8 bits)
    byteOffset1 += Uint8Array.BYTES_PER_ELEMENT;
    var extLength = getUint32(buffer, byteOffset);
    byteOffset1 += Uint32Array.BYTES_PER_ELEMENT;

    byteOffset1 += extLength;

    if (extId = 1) {

      exts.push( getNormals(buffer, byteOffset, extLength) );

    } else if (extId = 2) {

      exts.push( getWaterMask(buffer, byteOffset, extLength) );

    } else {
      console.warn(`ID ${extId} not recognised, only 1[= vertex normals] and 2[= water mask] are accepted`);
    }

    //getExtension();

  }

  return exts;

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

    let quantizedMeshComponents = {};

    let byteOffset = 0;

    let header = getHeader(buffer, byteOffset);
    byteOffset += 88;

    quantizedMeshComponents.header = header;

    var vertexCount = getUint32(buffer, byteOffset);
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


    var triangleCount = getUint32(buffer, byteOffset);
    byteOffset += Uint32Array.BYTES_PER_ELEMENT;

    //Immediately following the vertex data is the index data.
    //Indices specify how the vertices are linked together into triangles.
    //If tile has more than 65536 vertices, the tile uses the IndexData32 structure to encode indices.
    //Otherwise, it uses the IndexData16 structure.

    //To enforce proper byte alignment,
    //padding is added before the IndexData to ensure 2 byte alignment for IndexData16 and 4 byte alignment for IndexData32.

    var indices;

    if (vertexCount > 65536) {

      indices = getUint32Array(data, byteOffset, triangleCount * 3);
      byteOffset += triangleCount * 3 * 4;

    } else {

      indices = getUint16Array(data, byteOffset, triangleCount * 3);
      byteOffset += triangleCount * 3 * 2;

    }

    let indexArray = highwaterDecode(indices);

/////////////////////////////////

    let vertices = getVertices();

    let faces = getFaces();

    quantizedMeshComponents.geometry = getGeometry();

/////////////////////////////////

    // triangleCount * 3 = vertex count?
    var edgeIndicesArray = edgeIndicesDecode(buffer, byteOffset, triangleCount * 3);
    for (var i = 0, i < edgeIndicesArray, i++) {
      byteOffest += edgeIndicesArray[i].numBytes;
    }
    //byteOffset += edgeIndicesArray.west.length +  length

    quantizedMeshComponents.edges = returned edges objects

    //quantizedMeshComponents.extensions = [];

    quantizedMeshComponents.extensions = getExtension();

    //return new ThreeQuantizedMeshTile(header, uArray, vArray, heightArray, indexArray);

    //return Promise.resolve(quantizedMeshTile_bufferGeometry);

    return Promise.resolve(quantizedMeshComponents);

    //It should return a buffer geometry as part of an object containing vertexes, indices, edges, extensions

  }
};
