# JSON Referencing Test Suite

This repository contains a set of JSON objects that implementers of JSON referencing specifications can use to test their implementations.

It is meant to be language agnostic and should require only a JSON parser.
The conversion of the JSON objects into tests within a specific language and test framework of choice is left to be done by the implementer.

This suite is inspired by the [official JSON Schema Test Suite](https://github.com/json-schema-org/JSON-Schema-Test-Suite), where some of its tests originated.
Indeed JSON referencing is heavily influenced by JSON Schema, and it is only [recently](https://github.com/json-schema-org/referencing) that discussions have begun to formalize JSON referencing in a more cross-specification-amenable way.

## Structure of the Suite

The `tests` directory contains a set of folders corresponding to each specification which is tested by this suite.

Currently, this covers all modern JSON Schema specifications (notably, not yet OpenAPI specifications).
A `specifications.json` file is also included which maps each folder to a URL which identifies the specification (for JSON Schema these are known as "dialect ID"s).

Within each directory are tests corresponding to the particular specification.
Below is an example of such a test file, followed by a description of how to interpret the test.

```json
{
  "$schema": "../../test-schema.json",
  "registry": {
    "http://example.com/": {
      "definitions": {
        "foo": {
          "$id": "#foo",
          "foo": "bar"
        }
      }
    }
  },
  "tests": [
    {
      "base_uri": "http://example.com/",
      "ref": "#foo",
      "target": {
        "$id": "#foo",
        "foo": "bar"
      }
    }
  ]
}
```

Ignore the `$schema` property, it simply denotes that each test file satisfies a JSON Schema found at the given path.
The `registry` property contains a mapping between URIs and documents which are expected to be available for the duration of the tests.
The tests are found in the `tests` array, and each object within the array contain:

  * a `ref` key which is a reference to be resolved, along with an optional base URI to use to resolve a relative reference against

  * *either* of a `target` key which is the expected result of resolving the reference (taking into account the registry), *or* contain the key `error`, indicating that resolving the reference should produce some sort of error (because the reference is broken or somehow invalid)

  * an optional `then` key, which itself is a further test (recursively), and which is meant to be resolved *statefully* given the result of parent tests
