import array
import collections
import contextlib
import gc
import io
import pickle
import random
import signal
import struct
import sys
import tempfile
import tracemalloc
import unittest
from pathlib import Path
from unittest import mock


MODULE_DIR = Path(__file__).resolve().parents[1]
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import barcode_assembler as assembler_module
from barcode_assembler import HeadIndexQuery, ProjectAssembler


class _FakeBarcodeQuery:
    def __init__(self, values):
        self.values = values
        self.calls = collections.Counter()

    def get(self, read_int):
        self.calls[read_int] += 1
        return self.values.get(read_int, array.array("I"))


class _FakeHeadQuery:
    def __init__(self, hits):
        self.hits = hits
        self.calls = collections.Counter()

    def get(self, kmer):
        key = tuple(kmer)
        self.calls[key] += 1
        return self.hits.get(key, ())

    def cache_stats(self):
        return {
            "hits": 0,
            "misses": sum(self.calls.values()),
            "evictions": 0,
            "puts": 0,
            "rejections": 0,
            "entries": 0,
            "current_bytes": 0,
            "max_bytes": 0,
        }


def _graph_assembler(graph, **kwargs):
    """Represent a graph with one-barcode reads and a one-barcode head index."""
    nodes = set(graph)
    for candidates in graph.values():
        nodes.update(candidates)
    barcode_query = _FakeBarcodeQuery(
        {node: array.array("I", [node]) for node in nodes}
    )
    head_query = _FakeHeadQuery(
        {(node,): tuple(candidates) for node, candidates in graph.items()}
    )
    instance = ProjectAssembler(
        barcode_query=barcode_query,
        head_query=head_query,
        k_mer_size=1,
        end_fraction=1.0,
        end_cap=1,
        query_stride=1,
        **kwargs,
    )
    return instance, barcode_query, head_query


def _legacy_paths(graph, anchor_ids, max_depth, top_k, expansion_cap):
    """Literal reference for the path semantics in the pre-compact assembler."""
    anchor_ids = list(anchor_ids)
    targets = set(anchor_ids)
    results = {}
    for source_id in anchor_ids:
        queue = collections.deque(
            [(source_id, [source_id], frozenset([source_id]))]
        )
        seen_paths = {(source_id,)}
        source_paths = []
        expanded = 0
        capped = False
        while queue and not capped:
            _, path, path_set = queue.popleft()
            if len(path) > max_depth:
                continue
            for next_id in set(graph.get(path[-1], ())):
                if next_id in path_set:
                    continue
                next_path = path + [next_id]
                path_key = tuple(next_path)
                if path_key in seen_paths:
                    continue
                seen_paths.add(path_key)
                if next_id in targets and next_id != source_id:
                    source_paths.append(next_path)
                if len(next_path) < max_depth:
                    queue.append(
                        (next_id, next_path, path_set | {next_id})
                    )
                expanded += 1
                if expansion_cap is not None and expanded >= expansion_cap:
                    capped = True
                    break
        source_paths.sort(key=lambda path: (-len(path), tuple(path)))
        results[source_id] = source_paths[:top_k]
    return results


class CompactSearchTests(unittest.TestCase):
    def setUp(self):
        self.log_patcher = mock.patch.object(assembler_module, "log")
        self.log_patcher.start()
        self.addCleanup(self.log_patcher.stop)

    def assert_matches_legacy(
        self, graph, anchors, max_depth=5, top_k=10, cap=1000
    ):
        instance, _, _ = _graph_assembler(graph)
        actual = instance.find_top_paths_per_anchor(
            anchors,
            max_depth=max_depth,
            top_k_per_anchor=top_k,
            max_expansions_per_source=cap,
        )
        expected = _legacy_paths(graph, anchors, max_depth, top_k, cap)
        self.assertEqual(actual, expected)

    def test_branching_cycles_target_revisits_and_depth_match_legacy(self):
        graph = {
            0: (1, 2, 5, 6),
            1: (0, 3, 5),
            2: (3, 7),
            3: (5, 8, 9),
            5: (6,),
            6: (8,),
            7: (9,),
            8: (0,),
        }
        anchors = [0, 5, 6, 7, 8, 9]
        self.assert_matches_legacy(graph, anchors, max_depth=4)

    def test_random_small_graphs_match_legacy(self):
        generator = random.Random(3419)
        for case in range(30):
            graph = {
                node: tuple(
                    candidate
                    for candidate in range(8)
                    if generator.random() < 0.30
                )
                for node in range(8)
            }
            anchors = [0, 2, 5, 7]
            max_depth = generator.randrange(2, 6)
            top_k = generator.randrange(1, 8)
            cap = generator.randrange(1, 60)
            with self.subTest(case=case, depth=max_depth, cap=cap):
                self.assert_matches_legacy(
                    graph, anchors, max_depth, top_k, cap
                )

    def test_cap_is_counted_at_the_legacy_boundary(self):
        graph = {0: (1, 2), 1: (3,), 2: (3,)}
        self.assert_matches_legacy(
            graph, [0, 1, 2, 3], max_depth=4, top_k=10, cap=1
        )

    def test_cycles_do_not_consume_cap_and_target_is_traversed(self):
        # 0 -> 1 is retained because 1 is a target.  The back-edge 1 -> 0
        # must not consume the final expansion, so target 2 is retained too.
        graph = {0: (1,), 1: (0, 2)}
        instance, _, _ = _graph_assembler(graph)
        paths = instance.find_top_paths_per_anchor(
            [0, 1, 2],
            max_depth=3,
            top_k_per_anchor=10,
            max_expansions_per_source=2,
        )[0]
        self.assertEqual(paths, [[0, 1, 2], [0, 1]])

    def test_target_at_max_depth_and_legacy_depth_one_edge_are_retained(self):
        graph = {0: (1,), 1: (2,)}
        self.assert_matches_legacy(
            graph, [0, 2], max_depth=3, top_k=10, cap=10
        )
        self.assert_matches_legacy(
            graph, [0, 1], max_depth=1, top_k=10, cap=10
        )

    def test_ranking_is_longest_then_lexicographic_and_only_top_ten_return(self):
        targets = tuple(range(100, 130))
        graph = {0: (1,) + targets, 1: (2,), 2: targets}
        anchors = [0] + list(targets)
        instance, _, _ = _graph_assembler(graph)
        paths = instance.find_top_paths_per_anchor(
            anchors,
            max_depth=4,
            top_k_per_anchor=10,
            max_expansions_per_source=10000,
        )[0]
        self.assertEqual(
            paths,
            [[0, 1, 2, target] for target in targets[:10]],
        )

    def test_candidate_cache_preserves_deduplicated_iteration_order(self):
        barcodes = array.array("I", [1, 2])
        barcode_query = _FakeBarcodeQuery({7: barcodes})
        head_query = _FakeHeadQuery(
            {(1,): (9, 3, 9, 5), (2,): (5, 8, 3)}
        )
        instance = ProjectAssembler(
            barcode_query,
            head_query,
            k_mer_size=1,
            end_fraction=1.0,
            end_cap=2,
            query_stride=1,
            candidate_cache_mib=1,
        )
        legacy_candidates = set()
        legacy_candidates.update((9, 3, 9, 5))
        legacy_candidates.update((5, 8, 3))
        expected = list(legacy_candidates)
        first = instance._candidate_ids_for_read(7)
        second = instance._candidate_ids_for_read(7)
        self.assertEqual(list(first), expected)
        self.assertEqual(list(second), expected)
        self.assertEqual(barcode_query.calls[7], 1)
        self.assertEqual(sum(head_query.calls.values()), 2)

    def test_250k_expansion_frontier_uses_less_than_64_mib_traced_memory(self):
        class _WideHeadQuery:
            def get(self, kmer):
                if tuple(kmer) == (0,):
                    return range(1, 250001)
                return ()

            def cache_stats(self):
                return {
                    "hits": 0,
                    "misses": 1,
                    "evictions": 0,
                    "puts": 0,
                    "rejections": 0,
                    "entries": 0,
                    "current_bytes": 0,
                    "max_bytes": 0,
                }

        instance = ProjectAssembler(
            _FakeBarcodeQuery({0: array.array("I", [0])}),
            _WideHeadQuery(),
            k_mer_size=1,
            end_fraction=1.0,
            end_cap=1,
            query_stride=1,
            candidate_cache_mib=0,
            progress_every=0,
        )
        gc.collect()
        tracemalloc.start()
        try:
            result = instance.find_top_paths_per_anchor(
                [0],
                max_depth=3,
                top_k_per_anchor=10,
                max_expansions_per_source=250000,
            )
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        self.assertEqual(result, {0: []})
        self.assertLess(peak, 64 * 1024 * 1024)


class CacheBudgetTests(unittest.TestCase):
    def test_byte_lru_evicts_to_budget_and_rejects_oversized_entries(self):
        cache = assembler_module.ByteLRUCache(10)
        self.assertTrue(cache.put("a", "1234", 4))
        self.assertTrue(cache.put("b", "567890", 6))
        self.assertEqual(cache.get("a"), "1234")
        self.assertTrue(cache.put("c", "xx", 2))
        self.assertIsNone(cache.get("b"))
        self.assertLessEqual(cache.current_bytes, cache.max_bytes)
        self.assertFalse(cache.put("large", "x" * 11, 11))
        self.assertIsNone(cache.get("large"))
        stats = cache.stats()
        self.assertGreaterEqual(stats["hits"], 1)
        self.assertGreaterEqual(stats["misses"], 2)
        self.assertGreaterEqual(stats["evictions"], 1)

    def test_empty_candidate_entries_still_consume_budget(self):
        barcode_query = _FakeBarcodeQuery({})
        head_query = _FakeHeadQuery({})
        instance = ProjectAssembler(
            barcode_query,
            head_query,
            candidate_cache_mib=1,
            progress_every=0,
        )
        for read_int in range(10000):
            self.assertEqual(list(instance._candidate_ids_for_read(read_int)), [])
        stats = instance.candidate_cache_stats()
        self.assertLessEqual(stats["current_bytes"], stats["max_bytes"])
        self.assertLess(stats["entries"], 10000)
        self.assertGreater(stats["evictions"], 0)


class HeadIndexTests(unittest.TestCase):
    def _write_index(self, directory):
        keys_and_postings = [
            (tuple(range(11)), [9, 2, 9, 0x0102030405060708]),
            (tuple(range(11, 22)), [4]),
        ]
        idx = bytearray()
        data = bytearray()
        for key, postings in keys_and_postings:
            idx.extend(struct.pack(">11IQI", *key, len(data), len(postings)))
            data.extend(struct.pack(f">{len(postings)}Q", *postings))
        idx_path = directory / "head_index.idx"
        data_path = directory / "head_index.data"
        idx_path.write_bytes(idx)
        data_path.write_bytes(data)
        return idx_path, data_path, keys_and_postings

    def test_packed_key_lookup_posting_endianness_and_order(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            idx_path, data_path, entries = self._write_index(directory)
            query = HeadIndexQuery(idx_path, data_path)
            try:
                key, expected = entries[0]
                first = query.get(key)
                self.assertEqual(list(first), expected)
                first[0] = 12345
                self.assertEqual(list(query.get(key)), expected)
                self.assertEqual(list(query.get((99,) * 11)), [])
            finally:
                query.close()

    def test_metadata_cache_respects_byte_budget_and_evicts(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            idx_path, data_path, entries = self._write_index(directory)
            query = HeadIndexQuery(idx_path, data_path, cache_bytes=500)
            try:
                for key, _ in entries:
                    query.get(key)
                info = query.cache_stats()
                self.assertLessEqual(info["current_bytes"], 500)
                self.assertGreaterEqual(info["evictions"], 1)
            finally:
                query.close()


class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)
        self.input_file = self.directory / "anchors.pkl"
        self.input_file.write_bytes(b"anchors")
        self.anchors = [7, 3, 11]
        self.search_config = {
            "max_depth": 15,
            "top_k_per_anchor": 10,
            "max_expansions_per_source": 10000000,
            "k_mer_size": 11,
            "end_fraction": 0.25,
            "end_cap": 4000,
            "query_stride": 1,
        }
        self.metadata = assembler_module.build_checkpoint_metadata(
            self.anchors, [self.input_file], self.search_config
        )

    def tearDown(self):
        self.temporary.cleanup()

    def test_metadata_records_ordered_anchors_inputs_config_and_python(self):
        metadata = self.metadata
        self.assertEqual(
            metadata["schema_version"],
            assembler_module.CHECKPOINT_SCHEMA_VERSION,
        )
        self.assertEqual(metadata["python_version"], tuple(sys.version_info[:2]))
        self.assertEqual(metadata["anchor_count"], 3)
        self.assertEqual(metadata["search_config"], self.search_config)
        signature = next(iter(metadata["input_signatures"].values()))
        self.assertEqual(signature["size"], len(b"anchors"))
        self.assertIn("mtime_ns", signature)
        reordered = assembler_module.build_checkpoint_metadata(
            list(reversed(self.anchors)), [self.input_file], self.search_config
        )
        self.assertNotEqual(
            metadata["anchor_fingerprint"], reordered["anchor_fingerprint"]
        )

    def test_checkpoint_round_trip_and_exact_prefix_validation(self):
        checkpoint = self.directory / "bridges.checkpoint.pkl"
        paths = {7: [[7, 2, 3]], 3: []}
        assembler_module.write_checkpoint(
            checkpoint, self.metadata, [7, 3], paths
        )
        loaded = assembler_module.load_checkpoint(checkpoint, self.metadata)
        self.assertEqual(loaded["completed_anchor_ids"], [7, 3])
        self.assertEqual(loaded["paths_by_anchor"], paths)

        with self.assertRaises(assembler_module.CheckpointError):
            assembler_module.write_checkpoint(
                checkpoint, self.metadata, [7, 11], {7: [], 11: []}
            )

    def test_mismatched_or_corrupt_checkpoint_is_rejected_with_restart_hint(self):
        checkpoint = self.directory / "bridges.checkpoint.pkl"
        assembler_module.write_checkpoint(
            checkpoint, self.metadata, [7], {7: []}
        )
        different = dict(self.metadata)
        different["search_config"] = dict(self.search_config, max_depth=14)
        with self.assertRaisesRegex(
            assembler_module.CheckpointError, "--restart"
        ):
            assembler_module.load_checkpoint(checkpoint, different)

        checkpoint.write_bytes(b"not a pickle")
        with self.assertRaisesRegex(
            assembler_module.CheckpointError, "--restart"
        ):
            assembler_module.load_checkpoint(checkpoint, self.metadata)

    def test_atomic_pickle_preserves_old_file_when_serialization_fails(self):
        destination = self.directory / "bridges.pkl"
        assembler_module.atomic_write_pickle({"old": True}, destination)
        with mock.patch.object(
            assembler_module.pkl,
            "dump",
            side_effect=RuntimeError("synthetic failure"),
        ):
            with self.assertRaises(RuntimeError):
                assembler_module.atomic_write_pickle({"new": True}, destination)
        with destination.open("rb") as stream:
            self.assertEqual(pickle.load(stream), {"old": True})
        temporary_files = list(
            destination.parent.glob(f".{destination.name}.*.tmp")
        )
        self.assertEqual(temporary_files, [])

    def test_final_bridge_pickle_schema_is_unchanged(self):
        destination = self.directory / "bridges.pkl"
        bridges = {
            "read-a": [["read-a", "read-x", "read-b"]],
            "read-b": [],
        }
        assembler_module.atomic_write_pickle(bridges, destination)
        with destination.open("rb") as stream:
            loaded = pickle.load(stream)
        self.assertEqual(loaded, bridges)
        self.assertTrue(all(isinstance(key, str) for key in loaded))
        self.assertTrue(
            all(
                isinstance(read_id, str)
                for paths in loaded.values()
                for path in paths
                for read_id in path
            )
        )


class MiniaturePipelineTests(unittest.TestCase):
    @staticmethod
    def _write_pickle(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("wb") as stream:
            pickle.dump(value, stream, protocol=4)

    def _build_project(self, base_dir):
        self._write_pickle(
            base_dir / "rid_maps" / "rid_to_int.pkl",
            {"read-a": 0, "read-b": 1},
        )
        self._write_pickle(
            base_dir / "rid_maps" / "int_to_rid.pkl",
            {0: "read-a", 1: "read-b"},
        )
        self._write_pickle(
            base_dir / "output" / "anchors.pkl",
            {"read-a": {}, "read-b": {}},
        )

        store_dir = base_dir / "barcode_21mers"
        store_dir.mkdir(parents=True, exist_ok=True)
        first = array.array("I", range(44))
        second = array.array("I", range(100, 144))
        first_raw = first.tobytes()
        second_raw = second.tobytes()
        (store_dir / "read_unikmer_map.data").write_bytes(
            first_raw + second_raw
        )
        (store_dir / "read_unikmer_map.idx").write_bytes(
            struct.pack(">QII", 0, len(first), 1000)
            + struct.pack(">QII", len(first_raw), len(second), 1000)
        )

        head_dir = store_dir / "head_index"
        head_dir.mkdir()
        key = tuple(first[-11:])
        (head_dir / "head_index.idx").write_bytes(
            struct.pack(">11IQI", *key, 0, 1)
        )
        (head_dir / "head_index.data").write_bytes(struct.pack(">Q", 1))

    def test_end_to_end_binary_fixture_and_completed_checkpoint_resume(self):
        with tempfile.TemporaryDirectory() as temporary:
            base_dir = Path(temporary)
            self._build_project(base_dir)
            argv = [
                "--base-dir",
                str(base_dir),
                "--max-expansions-per-source",
                "100",
                "--candidate-cache-mib",
                "1",
                "--head-cache-mib",
                "1",
                "--progress-every",
                "0",
            ]
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                exit_status = assembler_module.main(argv)
            self.assertIn(exit_status, (None, 0))
            bridges_path = base_dir / "bridges_v1.1" / "bridges.pkl"
            with bridges_path.open("rb") as stream:
                first_run = pickle.load(stream)
            self.assertEqual(
                first_run,
                {
                    "read-a": [["read-a", "read-b"]],
                    "read-b": [],
                },
            )

            bridges_path.unlink()
            resumed_output = io.StringIO()
            with contextlib.redirect_stdout(resumed_output):
                exit_status = assembler_module.main(argv)
            self.assertIn(exit_status, (None, 0))
            with bridges_path.open("rb") as stream:
                resumed = pickle.load(stream)
            self.assertEqual(resumed, first_run)
            self.assertIn("RESUMED", resumed_output.getvalue())

    def test_interrupted_run_resumes_after_last_atomic_anchor_checkpoint(self):
        with tempfile.TemporaryDirectory() as temporary:
            base_dir = Path(temporary)
            self._build_project(base_dir)
            argv = [
                "--base-dir",
                str(base_dir),
                "--max-expansions-per-source",
                "100",
                "--candidate-cache-mib",
                "1",
                "--head-cache-mib",
                "1",
                "--progress-every",
                "0",
            ]
            original_search = ProjectAssembler.find_top_paths_per_anchor

            def interrupt_after_checkpoint(instance, *args, **kwargs):
                checkpoint_callback = kwargs["on_anchor_complete"]

                def stopping_callback(*callback_args):
                    checkpoint_callback(*callback_args)
                    raise assembler_module.TerminationRequested(signal.SIGTERM)

                kwargs["on_anchor_complete"] = stopping_callback
                return original_search(instance, *args, **kwargs)

            interrupted_output = io.StringIO()
            with mock.patch.object(
                ProjectAssembler,
                "find_top_paths_per_anchor",
                new=interrupt_after_checkpoint,
            ), contextlib.redirect_stdout(interrupted_output):
                exit_status = assembler_module.main(argv)
            self.assertEqual(exit_status, 128 + signal.SIGTERM)
            checkpoint_path = (
                base_dir
                / "bridges_v1.1"
                / "barcode_assembler.checkpoint.pkl"
            )
            with checkpoint_path.open("rb") as stream:
                checkpoint = pickle.load(stream)
            self.assertEqual(checkpoint["completed_anchor_ids"], [0])
            self.assertEqual(checkpoint["paths_by_anchor"], {0: [[0, 1]]})

            resumed_output = io.StringIO()
            with contextlib.redirect_stdout(resumed_output):
                exit_status = assembler_module.main(argv)
            self.assertEqual(exit_status, 0)
            with (base_dir / "bridges_v1.1" / "bridges.pkl").open(
                "rb"
            ) as stream:
                bridges = pickle.load(stream)
            self.assertEqual(
                bridges,
                {
                    "read-a": [["read-a", "read-b"]],
                    "read-b": [],
                },
            )
            self.assertIn("ANCHOR RESUMED 1/2 read-a", resumed_output.getvalue())


if __name__ == "__main__":
    unittest.main()
